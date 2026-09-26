#!/usr/bin/env python3
"""Verify and analyze archived sparse-candidate counts; never run design solves.

The standard-library analyzer consumes a compact bundle of the exact original
CSV/manifest bytes and frozen MMseqs2 cluster assignments. Optional plots require
matplotlib==3.10.6. All estimates describe the length-balanced sampled panel.
"""
import argparse
from collections import Counter, defaultdict
import csv
import fcntl
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import statistics
import sys
import tempfile

SCHEMA = "publication-density-v1"
BOOTSTRAP_SCHEMA = "paired-whole-cluster-stratum-calibration-v1"
HOST_NAMES = {"human":"Human", "yeast":"Yeast", "mouse":"Mouse", "zebrafish":"Zebrafish", "drosophila":"Drosophila", "celegans":"C. elegans", "arabidopsis":"Arabidopsis", "ecoli-o157-edl933":"E. coli O157:H7", "bsubtilis-species":"B. subtilis", "paeruginosa-pao1":"P. aeruginosa"}
HOST_ORDER = list(HOST_NAMES)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def file_bytes(path):
    raw = Path(path).read_bytes()
    return gzip.decompress(raw) if str(path).endswith(".gz") else raw


def read_csv(path):
    return list(csv.DictReader(io.StringIO(file_bytes(path).decode())))


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def save_csv(path, rows):
    if not rows:
        raise ValueError("Refusing to write an empty table")
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    raw = output.getvalue().encode()
    Path(path).write_bytes(gzip.compress(raw, mtime=0) if str(path).endswith(".gz") else raw)


def quantile(values, p):
    values = sorted(values)
    x = (len(values)-1)*p; lo = int(x); hi = min(lo+1, len(values)-1)
    return values[lo] + (x-lo)*(values[hi]-values[lo])


def key(row):
    return row["accession"], row["sequence_sha256"]


def prepare(root, out):
    campaign = root / "research/results/sparse-balanced-codon-panel-2026-07-16"
    inputs = out / "inputs"; (inputs / "merged").mkdir(parents=True, exist_ok=True)
    sources = {"plan-final.json.gz":campaign/"plan-final.json", "merged/manifest.json":campaign/"merged/manifest.json", "panel-2000.csv.gz":root/"research/data/uniprot-sprot-2026_02-sparse-balanced-2000.manifest.csv", "panel-400.csv.gz":root/"research/data/uniprot-sprot-2026_02-sparse-balanced-400.manifest.csv"}
    manifest = json.loads((campaign/"merged/manifest.json").read_text())
    for product in manifest["products"]:
        name = Path(product["path"]).name
        if name != product["objective_id"]+".csv":
            raise ValueError("Unexpected archived product name")
        sources["merged/"+name+".gz"] = campaign/"merged"/name
    records = []
    for rel, source in sources.items():
        raw = source.read_bytes()
        target = inputs / rel
        target.write_bytes(gzip.compress(raw, mtime=0) if rel.endswith(".gz") else raw)
        records.append({"path":rel, "uncompressed_sha256":sha(raw), "stored_sha256":sha(target.read_bytes()), "original_path":str(source.relative_to(root)), "uncompressed_bytes":len(raw)})
    save_json(inputs/"bundle.json", {"schema":SCHEMA, "files":records, "compression":"gzip with mtime=0; decompressed bytes equal original archived files", "scope":"Existing candidate counts and provenance only; no new sequence designs."})


def verified_inputs(out):
    inputs = out/"inputs"; bundle = json.loads((inputs/"bundle.json").read_text())
    for item in bundle["files"]:
        rel = Path(item["path"])
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError("Unsafe bundle path")
        p = inputs/rel
        if sha(p.read_bytes()) != item["stored_sha256"] or sha(file_bytes(p)) != item["uncompressed_sha256"]:
            raise ValueError("Archived input hash drift: "+str(rel))
    plan = json.loads(file_bytes(inputs/"plan-final.json.gz"))
    original_plan = dict(plan); original_plan.pop("plan_sha256")
    canonical = (json.dumps(original_plan, sort_keys=True, separators=(",", ":"), ensure_ascii=False)+"\n").encode()
    if sha(canonical) != plan["plan_sha256"]:
        raise ValueError("Archived plan's canonical hash is invalid")
    manifest = json.loads((inputs/"merged/manifest.json").read_text())
    if manifest["plan_sha256"] != plan["plan_sha256"]:
        raise ValueError("Merged products refer to a different plan")
    panel_raw = file_bytes(inputs/"panel-2000.csv.gz")
    if sha(panel_raw) != plan["inputs"]["manifest_csv"]["sha256"]:
        raise ValueError("Panel manifest differs from the planned source")
    panel_rows = read_csv(inputs/"panel-2000.csv.gz")
    panel = {key(r):r for r in panel_rows}
    subset_rows = read_csv(inputs/"panel-400.csv.gz"); subset = {key(r) for r in subset_rows}
    if len(panel_rows) != 2000 or len(panel) != 2000 or len(subset_rows) != 400 or len(subset) != 400 or not subset <= panel.keys():
        raise ValueError("Panel/subset membership or uniqueness failure")
    tasks = defaultdict(set)
    for task in plan["tasks"]:
        k = key(task)
        if k in tasks[task["objective_id"]]:
            raise ValueError("Duplicate planned task")
        tasks[task["objective_id"]].add(k)
        if int(task["aa_length"]) != int(panel[k]["aa_length"]) or task["length_bin"] != panel[k]["length_bin"]:
            raise ValueError("Planned task identity differs from panel")
    products, checks = {}, []
    for product in manifest["products"]:
        objective = product["objective_id"]
        if not re.fullmatch("[a-z0-9-]+", objective) or objective in products:
            raise ValueError("Duplicate or unsafe objective ID")
        path = inputs/"merged"/(objective+".csv.gz")
        if sha(file_bytes(path)) != product["sha256"]:
            raise ValueError("Merged CSV differs from its archived hash")
        rows = read_csv(path); rowkeys = {key(r) for r in rows}
        expected = panel.keys() if objective.startswith("human-") else subset
        if len(rows) != product["rows"] or len(rowkeys) != len(rows) or rowkeys != tasks[objective] or rowkeys != expected:
            raise ValueError("Merged objective membership does not equal the plan")
        for r in rows:
            if r["status"] != "ok" or r["error"] or r["solver_sha256"] != plan["inputs"]["solver"]["sha256"] or r["codon_table_sha256"] != product["codon_table_sha256"] or r["lambda"] != product["lambda"] or int(r["threads"]) != 8:
                raise ValueError("Mixed or failed archived measurement")
            n,z,d = (int(r[k]) for k in ("rna_nt","candidates","direct"))
            if n != 3*int(r["aa_length"]) or n != 3*int(panel[key(r)]["aa_length"]) or not 0 < z <= d:
                raise ValueError("Invalid dimensions/counts")
            if not math.isclose(z/n**2,float(r["Z_over_n2"]),rel_tol=5e-9,abs_tol=5e-12) or not math.isclose(z/d,float(r["candidate_over_direct"]),rel_tol=5e-9,abs_tol=5e-12):
                raise ValueError("Recorded candidate ratio does not reproduce")
            if int(r["storage_bytes"]) != 24*int(r["lattice_nodes"])+16*int(r["capacity"]):
                raise ValueError("Candidate-allocation arithmetic failure")
        products[objective] = rows
        checks.append({"objective":objective,"rows":len(rows),"sha256":product["sha256"],"codon_table_sha256":product["codon_table_sha256"]})
    if sum(map(len,products.values())) != 7600 or set(products) != set(tasks):
        raise ValueError("Incomplete campaign")
    cluster_dir = out/"clustering"
    cluster_prov = json.loads((cluster_dir/"run-provenance.json").read_text())
    if sha((cluster_dir/"clusters.tsv").read_bytes()) != cluster_prov["clusters_tsv_sha256"]:
        raise ValueError("Cluster assignments differ from their provenance")
    index = read_csv(cluster_dir/"sequence-index.csv")
    if {key(r) for r in index} != set(panel) or len(index) != 2000:
        raise ValueError("Clustering used a different panel")
    cluster_input = json.loads((cluster_dir/"input-provenance.json").read_text())
    if cluster_input["source_zip_sha256"] != plan["inputs"]["source"]["sha256"] or sha((cluster_dir/"sequence-index.csv").read_bytes()) != cluster_input["sequence_index_sha256"] or sha(file_bytes(cluster_dir/"panel.fasta.gz")) != cluster_input["fasta_sha256"]:
        raise ValueError("Clustering sequence source hash drift")
    clusters = {}
    for line in (cluster_dir/"clusters.tsv").read_text().splitlines():
        rep, member = line.split("\t")
        if member in clusters:
            raise ValueError("Duplicate cluster membership")
        clusters[member] = rep
    accessions = {r["accession"] for r in panel_rows}
    if set(clusters) != accessions or not set(clusters.values()) <= accessions or any(clusters[rep] != rep for rep in clusters.values()):
        raise ValueError("Clusters do not partition the panel")
    return products, panel, subset, clusters, {"checks":checks,"solver_sha256":plan["inputs"]["solver"]["sha256"],"plan_sha256":plan["plan_sha256"],"bundle_sha256":sha((inputs/"bundle.json").read_bytes()),"cluster_provenance":cluster_prov}


def fit(rows, weights=None):
    weights = weights or [1.0]*len(rows)
    xs = [math.log(int(r["rna_nt"])) for r in rows]; ys = [math.log(int(r["candidates"])) for r in rows]
    sw = math.fsum(weights); mx = math.fsum(w*x for w,x in zip(weights,xs))/sw; my = math.fsum(w*y for w,y in zip(weights,ys))/sw
    xx = math.fsum(w*(x-mx)**2 for w,x in zip(weights,xs)); xy = math.fsum(w*(x-mx)*(y-my) for w,x,y in zip(weights,xs,ys)); yy = math.fsum(w*(y-my)**2 for w,y in zip(weights,ys))
    beta = xy/xx; intercept = my-beta*mx
    sse = math.fsum(w*(y-intercept-beta*x)**2 for w,x,y in zip(weights,xs,ys))
    return {"n":len(rows),"weight_sum":sw,"beta":beta,"intercept_ln":intercept,"coefficient":math.exp(intercept),"r_squared":1-sse/yy,"sse_log":sse,"rmse_log":math.sqrt(sse/sw)}


def linear_solve(matrix, rhs):
    a = [list(row)+[value] for row,value in zip(matrix,rhs)]
    for i in range(len(a)):
        j = max(range(i,len(a)),key=lambda j:abs(a[j][i])); a[i],a[j]=a[j],a[i]
        pivot=a[i][i]
        if abs(pivot)<1e-12: raise ValueError("Singular fit")
        a[i]=[v/pivot for v in a[i]]
        for j in range(len(a)):
            if i!=j:
                factor=a[j][i];a[j]=[u-factor*v for u,v in zip(a[j],a[i])]
    return [row[-1] for row in a]


def hinge_fit(rows, knot_aa):
    knot=math.log(3*knot_aa)
    design = [[1,math.log(int(r["rna_nt"])),max(0,math.log(int(r["rna_nt"]))-knot)] for r in rows]
    ys = [math.log(int(r["candidates"])) for r in rows]
    normal = [[math.fsum(v[i]*v[j] for v in design) for j in range(3)] for i in range(3)]
    rhs = [math.fsum(v[i]*y for v,y in zip(design,ys)) for i in range(3)]
    a,b,c = linear_solve(normal,rhs)
    sse = math.fsum((y-a-b*v[1]-c*v[2])**2 for v,y in zip(design,ys)); plain=fit(rows)
    return {"knot_aa":knot_aa,"intercept_ln":a,"lower_slope":b,"upper_slope":b+c,"hinge_coefficient":c,"sse_log":sse,"sse_reduction_fraction":1-sse/plain["sse_log"],"r_squared":1-(1-plain["r_squared"])*sse/plain["sse_log"],"n":len(rows)}


def describe(rows):
    ratios = [int(r["candidates"])/int(r["direct"]) for r in rows]
    density = [int(r["candidates"])/int(r["rna_nt"])**2 for r in rows]
    allocation = [int(r["storage_bytes"])/2**20 for r in rows]
    result={"rows":len(rows),"fit":fit(rows)}
    for name,values in (("retention",ratios),("density",density),("candidate_allocation_mib",allocation)):
        result[name]={label:quantile(values,p) for label,p in (("min",0),("p05",.05),("median",.5),("p95",.95),("max",1))}
    return result


def moment_slope(s, column):
    return (s[column+1]-s[1]*s[column]/s[0])/(s[2]-s[1]*s[1]/s[0])


def bootstrap(products, panel, clusters, replicates, seed):
    human4 = {key(r):r for r in products["human-lambda4"]}
    groups=defaultdict(list); cluster_members=defaultdict(list)
    for r in products["human-lambda0"]:
        x=math.log(int(r["rna_nt"]));y0=math.log(int(r["candidates"]));y4=math.log(int(human4[key(r)]["candidates"]))
        vector=(1.,x,x*x,y0,x*y0,y4,x*y4)
        groups[panel[key(r)]["length_bin"]].append(vector)
        cluster_members[clusters[r["accession"]]].append((panel[key(r)]["length_bin"],vector))
    labels=sorted(groups); sizes=[len(groups[h]) for h in labels]; label_index={h:i for i,h in enumerate(labels)}
    cache=[]
    for rep in sorted(cluster_members):
        group=cluster_members[rep]; by_stratum=defaultdict(lambda:[0.]*7)
        for h,v in group:
            s=by_stratum[label_index[h]]
            for j in range(7):s[j]+=v[j]
        cache.append([(h,tuple(s),tuple(v/len(group) for v in s)) for h,s in sorted(by_stratum.items())])
    rng_within=random.Random(seed);rng_cluster=random.Random(seed+1)
    output=[]
    for iteration in range(replicates):
        sample=[v for h in labels for v in rng_within.choices(groups[h],k=len(groups[h]))]
        s=[math.fsum(v[j] for v in sample) for j in range(7)]
        record={"replicate":iteration+1,"within_lambda0":moment_slope(s,3),"within_lambda4":moment_slope(s,5)}
        strata=[[0.]*7 for _ in labels];weighted=[[0.]*7 for _ in labels]
        for cluster in rng_cluster.choices(cache,k=len(cache)):
            for h,plain,weight in cluster:
                a=strata[h];b=weighted[h]
                for j in range(7):a[j]+=plain[j];b[j]+=weight[j]
        if any(s[0]==0 for s in strata):
            raise ValueError("Whole-cluster bootstrap lost an entire length stratum; prespecify a revised procedure")
        calibrated=[math.fsum(s[j]*target/s[0] for s,target in zip(strata,sizes)) for j in range(7)]
        equal=[math.fsum(s[j] for s in weighted) for j in range(7)]
        weighted_calibrated=[math.fsum(s[j]*target/s[0] for s,target in zip(weighted,sizes)) for j in range(7)]
        for name,v in (("cluster",calibrated),("equal_cluster",equal),("calibrated_equal_cluster",weighted_calibrated)):
            record[name+"_lambda0"]=moment_slope(v,3);record[name+"_lambda4"]=moment_slope(v,5)
        output.append(record)
        if (iteration+1)%1000==0:print(f"Bootstrap {iteration+1}/{replicates}",flush=True)
    return output


def analyze(out, products, panel, subset, clusters, provenance, draws):
    human={};residuals=[];strata_rows=[];cutoff_rows=[];piecewise_rows=[];cluster_sizes=Counter(clusters.values())
    sizes=Counter(r["length_bin"] for r in panel.values())
    for lam in (0,4):
        objective=f"human-lambda{lam}";rows=products[objective]
        value=describe(rows);value["cutoffs"]={}
        for cutoff in (0,100,400,1000):
            f=fit([r for r in rows if int(r["aa_length"])>=cutoff]);value["cutoffs"][str(cutoff)]=f
            cutoff_rows.append(dict(lambda_user=lam,minimum_aa=cutoff,**f))
        value["disjoint_length_band_fits"]={}
        for lo,hi in ((10,99),(100,399),(400,999),(1000,3999)):
            selected=[r for r in rows if lo<=int(r["aa_length"])<=hi]
            f=fit(selected);value["disjoint_length_band_fits"][f"{lo}-{hi}"]=f
            piecewise_rows.append(dict(lambda_user=lam,minimum_aa=lo,maximum_aa=hi,**f))
        value["continuous_hinge_fits"]={str(k):hinge_fit(rows,k) for k in (100,400,1000)}
        weights=[1/cluster_sizes[clusters[r["accession"]]] for r in rows]
        totals=defaultdict(float)
        for r,w in zip(rows,weights):totals[panel[key(r)]["length_bin"]]+=w
        calibrated=[w*sizes[panel[key(r)]["length_bin"]]/totals[panel[key(r)]["length_bin"]] for r,w in zip(rows,weights)]
        value["equal_cluster_fit"]=fit(rows,weights);value["stratum_calibrated_equal_cluster_fit"]=fit(rows,calibrated)
        value["bootstrap_95_ci"]={name:[quantile([float(d[name+f"_lambda{lam}"]) for d in draws],p) for p in (.025,.975)] for name in ("within","cluster","equal_cluster","calibrated_equal_cluster")}
        for r in rows:
            x=math.log(int(r["rna_nt"]));y=math.log(int(r["candidates"]));f=value["fit"];h=value["continuous_hinge_fits"]["400"]
            residuals.append(dict(accession=r["accession"],lambda_user=lam,aa_length=int(r["aa_length"]),rna_nt=int(r["rna_nt"]),length_bin=panel[key(r)]["length_bin"],cluster=clusters[r["accession"]],candidates=int(r["candidates"]),retention=int(r["candidates"])/int(r["direct"]),density=int(r["candidates"])/int(r["rna_nt"])**2,log_n=x,log_z=y,linear_residual=y-f["intercept_ln"]-f["beta"]*x,hinge400_residual=y-h["intercept_ln"]-h["lower_slope"]*x-h["hinge_coefficient"]*max(0,x-math.log(1200))))
        for label in sorted(sizes):
            selected=[r for r in rows if panel[key(r)]["length_bin"]==label];s=describe(selected)
            res=[r["linear_residual"] for r in residuals if r["lambda_user"]==lam and r["length_bin"]==label]
            strata_rows.append(dict(lambda_user=lam,length_bin=label,n=len(selected),minimum_aa=min(int(r["aa_length"]) for r in selected),maximum_aa=max(int(r["aa_length"]) for r in selected),median_rna_nt=statistics.median(int(r["rna_nt"]) for r in selected),median_density=s["density"]["median"],p05_density=s["density"]["p05"],p95_density=s["density"]["p95"],median_retention=s["retention"]["median"],p05_retention=s["retention"]["p05"],p95_retention=s["retention"]["p95"],median_candidate_mib=s["candidate_allocation_mib"]["median"],p95_candidate_mib=s["candidate_allocation_mib"]["p95"],median_linear_residual=statistics.median(res),mean_linear_residual=statistics.mean(res)))
        human[objective]=value
    human0={key(r):r for r in products["human-lambda0"]};human4={key(r):r for r in products["human-lambda4"]}
    host={};host_rows=[]
    for name in HOST_ORDER:
        objective=name+"-lambda4"; rows=[r for r in products[objective] if key(r) in subset]
        value=describe(rows);value["median_paired_candidate_ratio_to_human"]=statistics.median(int(r["candidates"])/int(human4[key(r)]["candidates"]) for r in rows);value["codon_table_sha256"]=rows[0]["codon_table_sha256"]
        value["frequency_precision"]="rounded two-decimal synonymous frequencies" if name in ("human","yeast") else "full within-synonymous precision from archived Kazusa counts"
        host[objective]=value
        host_rows.append(dict(host=name,label=HOST_NAMES[name],n=len(rows),median_density=value["density"]["median"],p05_density=value["density"]["p05"],p95_density=value["density"]["p95"],median_retention=value["retention"]["median"],paired_ratio_to_human=value["median_paired_candidate_ratio_to_human"],median_candidate_mib=value["candidate_allocation_mib"]["median"],p95_candidate_mib=value["candidate_allocation_mib"]["p95"],beta=value["fit"]["beta"],r_squared=value["fit"]["r_squared"]))
    by_cluster=defaultdict(list)
    for k,row in panel.items():by_cluster[clusters[row["accession"]]].append(row)
    cluster_rows=[dict(representative=rep,size=len(members),minimum_aa=min(int(r["aa_length"]) for r in members),maximum_aa=max(int(r["aa_length"]) for r in members),strata=len({r["length_bin"] for r in members})) for rep,members in sorted(by_cluster.items())]
    result={"schema":SCHEMA,"source_script_sha256":sha(Path(__file__).read_bytes()),"scope":"Read-only analysis of archived candidate counts; no new optimization. These are unweighted length-balanced panel descriptions, not Swiss-Prot population estimates or complexity bounds.","total_rows":sum(map(len,products.values())),"human":human,"hosts":host,"paired_human":{"strictly_fewer_candidates_at_lambda4":sum(int(human4[k]["candidates"])<int(r["candidates"]) for k,r in human0.items()),"median_candidate_ratio":statistics.median(int(human4[k]["candidates"])/int(r["candidates"]) for k,r in human0.items()),"median_allocation_ratio":statistics.median(int(human4[k]["storage_bytes"])/int(r["storage_bytes"]) for k,r in human0.items()),"median_mfe_change_per_nt":statistics.median((float(human4[k]["mfe_kcal"])-float(r["mfe_kcal"]))/int(r["rna_nt"]) for k,r in human0.items()),"median_cai_lambda0":statistics.median(float(r["cai"]) for r in human0.values()),"median_cai_lambda4":statistics.median(float(r["cai"]) for r in human4.values())},"clusters":{"count":len(by_cluster),"singletons":sum(len(v)==1 for v in by_cluster.values()),"largest":max(map(len,by_cluster.values())),"size_distribution":dict(sorted(Counter(map(len,by_cluster.values())).items())),"cross_stratum_clusters":sum(len({r["length_bin"] for r in v})>1 for v in by_cluster.values()),"definition":"MMseqs2 operational sequence-similarity groups at 30% alignment identity and 80% bidirectional coverage, plus recorded heuristics/E-value; not asserted biological families."},"provenance":provenance,"bootstrap":{"replicates":len(draws),"within_stratum_seed":20260716,"whole_cluster_seed":20260717,"schema":BOOTSTRAP_SCHEMA,"calibration":"Whole-cluster resampling uses one multiplicity per cluster across both lambdas and all strata, then calibrates each resampled stratum to its original sequence count. Equal-cluster fits use weight 1/cluster-size; a separately reported variant also calibrates these weights within length strata."},"summed_instrumented_task_wall_hours":sum(float(r["wall_seconds"]) for rows in products.values() for r in rows)/3600,"limitations":["Archived CSVs lack emitted sequences/structures and do not establish global design optimality.","Candidate inventories and allocation are not actual split work or whole-process RSS.","Bootstrap intervals are conditional panel sensitivities; cluster assignments are heuristic and do not establish independent biological families.","Finite-range curvature and cutoff sensitivity preclude extrapolating a single power-law exponent.","Host fits use the same 400 proteins and report point fits only."]}
    tables=out/"tables";tables.mkdir(exist_ok=True)
    for name,rows in (("strata.csv",strata_rows),("cutoffs.csv",cutoff_rows),("disjoint-band-fits.csv",piecewise_rows),("hosts.csv",host_rows),("cluster-sizes.csv",cluster_rows),("residuals.csv.gz",residuals)):
        save_csv(tables/name,rows)
    save_json(out/"summary.json",result)
    write_tex(out,result,strata_rows,host_rows)
    write_results(out,result)
    return result


def write_tex(out,result,strata_rows,host_rows):
    macros={"DensityTaskCount":str(result["total_rows"]),"DensityProteinCount":str(result["human"]["human-lambda0"]["rows"]),"DensityCodonTableCount":str(len(result["hosts"])),"DensityHostSubsetCount":str(result["hosts"]["human-lambda4"]["rows"]),"DensityClusterCount":str(result["clusters"]["count"]),"DensitySingletonClusters":str(result["clusters"]["singletons"]),"DensityLargestCluster":str(result["clusters"]["largest"]),"DensityPairedRatio":f"{result['paired_human']['median_candidate_ratio']:.3f}","DensityPairedReducedCount":str(result["paired_human"]["strictly_fewer_candidates_at_lambda4"]),"DensityBootstrapReplicates":str(result["bootstrap"]["replicates"])}
    for lam,word in ((0,"Zero"),(4,"Four")):
        s=result["human"][f"human-lambda{lam}"];prefix="Density"+word
        macros.update({prefix+"Median":f"{s['density']['median']:.6f}",prefix+"RetentionPercent":f"{100*s['retention']['median']:.2f}",prefix+"MaxRetentionPercent":f"{100*s['retention']['max']:.2f}",prefix+"Beta":f"{s['fit']['beta']:.3f}",prefix+"RSquared":f"{s['fit']['r_squared']:.3f}",prefix+"ClusterBeta":f"{s['equal_cluster_fit']['beta']:.3f}",prefix+"CalibratedClusterBeta":f"{s['stratum_calibrated_equal_cluster_fit']['beta']:.3f}"})
        for name,label in (("within","Within"),("cluster","Cluster"),("equal_cluster","EqualCluster"),("calibrated_equal_cluster","CalibratedCluster")):
            lo,hi=s["bootstrap_95_ci"][name];macros[prefix+label+"CILow"]=f"{lo:.3f}";macros[prefix+label+"CIHigh"]=f"{hi:.3f}"
        for cutoff,label in ((100,"Hundred"),(400,"FourHundred"),(1000,"Thousand")):
            macros[prefix+label+"Beta"]=f"{s['cutoffs'][str(cutoff)]['beta']:.3f}"
        hinge=s["continuous_hinge_fits"]["400"]
        macros.update({prefix+"HingeLowerBeta":f"{hinge['lower_slope']:.3f}",prefix+"HingeUpperBeta":f"{hinge['upper_slope']:.3f}",prefix+"HingeReductionPercent":f"{100*hinge['sse_reduction_fraction']:.1f}"})
    betas=[h["fit"]["beta"] for h in result["hosts"].values()];medians=[h["density"]["median"] for h in result["hosts"].values()]
    macros.update(DensityHostBetaMinimum=f"{min(betas):.3f}",DensityHostBetaMaximum=f"{max(betas):.3f}",DensityHostMedianMinimum=f"{min(medians):.6f}",DensityHostMedianMaximum=f"{max(medians):.6f}")
    (out/"claims.tex").write_text("% Generated from verified archived inputs by publication_density.py\n"+"".join("\\newcommand{\\"+k+"}{"+v+"}\n" for k,v in macros.items()))
    lines=[r"\begin{tabular}{rrrr}",r"\toprule",r"Minimum aa & Proteins & $\beta$, $\lambda=0$ & $\beta$, $\lambda=4$ \\",r"\midrule"]
    for cutoff in (0,100,400,1000):
        a=result["human"]["human-lambda0"]["cutoffs"][str(cutoff)];b=result["human"]["human-lambda4"]["cutoffs"][str(cutoff)]
        lines.append(f"{10 if cutoff==0 else cutoff} & {a['n']:,} & {a['beta']:.3f} & {b['beta']:.3f} \\\\")
    lines += [r"\bottomrule",r"\end{tabular}"]
    (out/"tables/cutoffs.tex").write_text("\n".join(lines)+"\n")
    lines=[r"\begin{tabular}{lrrrr}",r"\toprule",r"Codon table & Median $Z/n^2$ & Ratio to human & $\beta$ & $R^2$ \\",r"\midrule"]
    for row in host_rows:lines.append(f"{row['label']} & {row['median_density']:.6f} & {row['paired_ratio_to_human']:.3f} & {row['beta']:.3f} & {row['r_squared']:.3f} \\\\")
    lines += [r"\bottomrule",r"\end{tabular}"]
    (out/"tables/hosts.tex").write_text("\n".join(lines)+"\n")
    lines=[r"\begin{tabular}{lrr}",r"\toprule",r"Fit / interval & $\lambda=0$ & $\lambda=4$ \\",r"\midrule"]
    for label,key,fitkey in (("Panel / within-stratum bootstrap","within","fit"),("Panel / whole-cluster bootstrap","cluster","fit"),("Equal-cluster weighting","equal_cluster","equal_cluster_fit"),("Equal-cluster, length calibrated","calibrated_equal_cluster","stratum_calibrated_equal_cluster_fit")):
        values=[]
        for lam in (0,4):
            s=result["human"][f"human-lambda{lam}"];lo,hi=s["bootstrap_95_ci"][key]
            values.append(f"{s[fitkey]['beta']:.3f} [{lo:.3f}, {hi:.3f}]")
        lines.append(label+" & "+" & ".join(values)+r" \\")
    lines += [r"\bottomrule",r"\end{tabular}"]
    (out/"tables/clusters.tex").write_text("\n".join(lines)+"\n")
    lines=[r"\begin{tabular}{lrrrrr}",r"\toprule",r"Length stratum (aa) & Proteins & $Z/n^2$, $\lambda=0$ & $Z/n^2$, $\lambda=4$ & Retention, $\lambda=0$ & Retention, $\lambda=4$ \\",r"\midrule"]
    for label in sorted({r["length_bin"] for r in strata_rows}):
        a=next(r for r in strata_rows if r["length_bin"]==label and r["lambda_user"]==0);b=next(r for r in strata_rows if r["length_bin"]==label and r["lambda_user"]==4)
        lo,hi=map(int,label.removeprefix("aa").split("-"))
        lines.append(f"{lo:,}--{hi:,} & {a['n']} & {a['median_density']:.6f} & {b['median_density']:.6f} & {100*a['median_retention']:.2f}\\% & {100*b['median_retention']:.2f}\\% \\\\")
    lines += [r"\bottomrule",r"\end{tabular}"]
    (out/"tables/strata.tex").write_text("\n".join(lines)+"\n")
    lines=[r"\begin{tabular}{lrr}",r"\toprule",r"Statistic & $\lambda=0$ & $\lambda=4$ \\",r"\midrule"]
    for label,name,q,scale,precision in (("Median $Z/n^2$","density","median",1,6),("P05 $Z/n^2$","density","p05",1,6),("P95 $Z/n^2$","density","p95",1,6),("Median retention (percent)","retention","median",100,2),("Maximum retention (percent)","retention","max",100,2),("Median candidate allocation (MiB)","candidate_allocation_mib","median",1,2),("P95 candidate allocation (MiB)","candidate_allocation_mib","p95",1,2),("Maximum candidate allocation (MiB)","candidate_allocation_mib","max",1,2)):
        values=[f"{scale*result['human'][f'human-lambda{lam}'][name][q]:.{precision}f}" for lam in (0,4)]
        lines.append(label+" & "+" & ".join(values)+r" \\")
    lines += [r"\bottomrule",r"\end{tabular}"]
    (out/"tables/human.tex").write_text("\n".join(lines)+"\n")


def write_results(out,s):
    a=s["human"]["human-lambda0"];b=s["human"]["human-lambda4"];c=s["clusters"]
    text=["# Generated archived-density findings", "", f"Verified all **{s['total_rows']:,}** archived records across {len(s['provenance']['checks'])} merged products, including original CSV hashes, the canonical campaign-plan hash, exact task/panel membership, sequence identities, provenance fields and candidate/storage arithmetic. This work performs no new design optimization.", "", f"The two human panels contain {a['rows']:,} proteins each. Median candidate retention is {100*a['retention']['median']:.4f}% at lambda 0 and {100*b['retention']['median']:.4f}% at lambda 4; maxima are {100*a['retention']['max']:.4f}% and {100*b['retention']['max']:.4f}%. Median Z/n² is {a['density']['median']:.8f} and {b['density']['median']:.8f}. All {s['paired_human']['strictly_fewer_candidates_at_lambda4']:,} paired proteins retain fewer candidates at lambda 4; the median ratio is {s['paired_human']['median_candidate_ratio']:.5f}.", "", "## Length sensitivity", "", "| Minimum included aa | Proteins | Slope, lambda 0 | Slope, lambda 4 |", "|---:|---:|---:|---:|"]
    for cutoff in (0,100,400,1000):
        aa=a["cutoffs"][str(cutoff)];bb=b["cutoffs"][str(cutoff)]
        text.append(f"| {10 if cutoff==0 else cutoff} | {aa['n']:,} | {aa['beta']:.6f} | {bb['beta']:.6f} |")
    text += ["", "The residual plot shows curvature. Continuous two-slope fits with a prespecified 400-aa knot give:", "", "| Lambda | Lower slope | Upper slope | Reduction in log-space SSE |", "|---:|---:|---:|---:|"]
    for lam,value in ((0,a),(4,b)):
        h=value["continuous_hinge_fits"]["400"];text.append(f"| {lam} | {h['lower_slope']:.6f} | {h['upper_slope']:.6f} | {100*h['sse_reduction_fraction']:.2f}% |")
    text += ["", "This is descriptive model sensitivity, not an asymptotic complexity estimate. Full tables also include disjoint length-band fits and fixed knots at 100 and 1000 aa.", "", "## Sequence-similarity sensitivity", "", f"Read-only MMseqs2 clustering yields {c['count']:,} operational similarity groups: {c['singletons']:,} singletons, maximum group size {c['largest']}, and {c['cross_stratum_clusters']} groups crossing length strata. The recorded settings require at least 30% alignment identity and 80% coverage of both sequences, with additional prefilter/E-value settings. These are heuristic algorithm-defined clusters, not verified biological families.", "", f"Each interval below uses {s['bootstrap']['replicates']:,} paired replicates. Whole-cluster resampling preserves one multiplicity per cluster across lambdas and strata, then calibrates the original stratum totals. Equal-cluster fits give each cluster total weight one; their length-calibrated variant preserves the original stratum masses.", "", "| Fit / bootstrap | Lambda 0: slope [95% interval] | Lambda 4: slope [95% interval] |", "|---|---:|---:|"]
    for label,name,f in (("Panel / within strata","within","fit"),("Panel / whole clusters","cluster","fit"),("Equal clusters","equal_cluster","equal_cluster_fit"),("Equal clusters + length strata","calibrated_equal_cluster","stratum_calibrated_equal_cluster_fit")):
        vals=[]
        for h in (a,b):
            lo,hi=h["bootstrap_95_ci"][name];vals.append(f"{h[f]['beta']:.6f} [{lo:.6f}, {hi:.6f}]")
        text.append("| "+label+" | "+" | ".join(vals)+" |")
    betas=[v["fit"]["beta"] for v in s["hosts"].values()]
    text += ["", "The length-range effect is larger than these cluster adjustments in this panel. Bootstrap intervals remain conditional sensitivity measures and cannot establish independent biological families or extrapolation validity.", "", "## Codon tables and artifacts", "", f"The {len(s['hosts'])} codon tables use the same {s['hosts']['human-lambda4']['rows']} proteins at lambda 4. Their point-fit slopes range from {min(betas):.6f} to {max(betas):.6f}. Human and yeast frequencies are rounded; the other archived tables retain full within-synonymous precision. Host fits report point estimates only.", "", "The figures `density-sensitivity.pdf` and `host-comparison.pdf`, generated `claims.tex`, and LaTeX/CSV tables derive from the verified data. The raw bootstrap draws, cluster assignments, compressed original inputs and exact source snapshots are retained for reproduction. Candidate allocation describes candidate containers, not whole-process RSS; archived instrumented wall time is not a controlled performance benchmark.", ""]
    (out/"RESULTS.md").write_text("\n".join(text))


def plots(out):
    os.environ.setdefault("MPLCONFIGDIR",tempfile.mkdtemp(prefix="sparsedesign-density-mpl-"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    s=json.loads((out/"summary.json").read_text());rows=read_csv(out/"tables/residuals.csv.gz");strata=read_csv(out/"tables/strata.csv");colors={0:"#17658c",4:"#b55931"}
    plt.rcParams.update({"font.size":8,"axes.labelsize":9,"axes.titlesize":10,"legend.fontsize":8,"pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
    fig,ax=plt.subplots(2,2,figsize=(7.1,5.3),layout="constrained")
    for lam in (0,4):
        rr=[r for r in rows if int(r["lambda_user"])==lam];ss=[r for r in strata if int(r["lambda_user"])==lam];human=s["human"][f"human-lambda{lam}"]
        x=[float(r["median_rna_nt"]) for r in ss]
        ax[0,0].plot(x,[100*float(r["median_retention"]) for r in ss],marker="o",markersize=3,color=colors[lam],label=rf"$\lambda={lam}$")
        ax[0,0].fill_between(x,[100*float(r["p05_retention"]) for r in ss],[100*float(r["p95_retention"]) for r in ss],color=colors[lam],alpha=.12)
        ax[0,1].plot([10,100,400,1000],[human["cutoffs"][str(k)]["beta"] for k in (0,100,400,1000)],marker="o",markersize=4,color=colors[lam],label=rf"$\lambda={lam}$")
        ax[1,0].scatter([int(r["rna_nt"]) for r in rr],[float(r["linear_residual"]) for r in rr],s=3,alpha=.14,color=colors[lam],rasterized=True)
        ax[1,0].plot(x,[float(r["median_linear_residual"]) for r in ss],color=colors[lam],marker="o",markersize=3,label=rf"$\lambda={lam}$ stratum median")
        for i,(key,fitkey) in enumerate((("within","fit"),("cluster","fit"),("equal_cluster","equal_cluster_fit"),("calibrated_equal_cluster","stratum_calibrated_equal_cluster_fit"))):
            point=human[fitkey]["beta"];lo,hi=human["bootstrap_95_ci"][key]
            ax[1,1].errorbar(point,i+(-.12 if lam==0 else .12),xerr=[[point-lo],[hi-point]],fmt="o",markersize=4,color=colors[lam],capsize=2)
    for a in (ax[0,0],ax[0,1],ax[1,0]):a.set_xscale("log");a.grid(axis="y",alpha=.18)
    ax[0,0].set_yscale("log");ax[0,0].set_ylabel("Candidate retention Z/direct (%)");ax[0,0].set_xlabel("RNA length n (nt)");ax[0,0].legend(frameon=False);ax[0,0].set_title("A  Length-stratum median and P05–P95",loc="left")
    ax[0,1].set_xticks([10,100,400,1000],labels=["10","100","400","1000"]);ax[0,1].set_xlabel("Minimum protein length included (aa)");ax[0,1].set_ylabel(r"Fitted log–log slope $\beta$");ax[0,1].set_title("B  Cutoff sensitivity",loc="left")
    ax[1,0].axhline(0,color="black",linewidth=.6,alpha=.6);ax[1,0].set_xlabel("RNA length n (nt)");ax[1,0].set_ylabel("Residual from global log–log fit");ax[1,0].set_title("C  Curvature beyond one fitted slope",loc="left")
    ax[1,1].set_yticks(range(4),labels=["Within strata","Whole clusters","Equal clusters","Equal clusters\n+ length strata"]);ax[1,1].invert_yaxis();ax[1,1].set_xlabel(r"Slope $\beta$ and percentile 95% interval");ax[1,1].set_title("D  Resampling / weighting sensitivity",loc="left");ax[1,1].grid(axis="x",alpha=.18)
    fig.savefig(out/"density-sensitivity.pdf");fig.savefig(out/"density-sensitivity.png",dpi=180);plt.close(fig)
    hosts=read_csv(out/"tables/hosts.csv");fig,ax=plt.subplots(1,2,figsize=(7.1,3.4),layout="constrained",sharey=True)
    for i,row in enumerate(hosts):
        med=float(row["median_density"]);lo=float(row["p05_density"]);hi=float(row["p95_density"])
        ax[0].errorbar(med,i,xerr=[[med-lo],[hi-med]],fmt="o",markersize=4,color="#17658c",capsize=2)
        ax[1].plot(float(row["beta"]),i,"o",markersize=4,color="#17658c")
    ax[0].set_xscale("log");ax[0].set_yticks(range(len(hosts)),labels=[r["label"] for r in hosts]);ax[0].invert_yaxis();ax[0].set_xlabel(r"Median $Z/n^2$, with P05–P95");ax[0].set_title("A  Same 400-protein panel",loc="left")
    ax[1].axvline(float(hosts[0]["beta"]),color="grey",linestyle="--",linewidth=.8);ax[1].set_xlabel(r"Point-fit slope $\beta$ at $\lambda=4$");ax[1].set_title("B  Codon-table comparison",loc="left")
    for a in ax:a.grid(axis="x",alpha=.18)
    fig.savefig(out/"host-comparison.pdf");fig.savefig(out/"host-comparison.png",dpi=180);plt.close(fig)
    save_json(out/"plot-provenance.json",{"matplotlib":matplotlib.__version__,"script_sha256":sha(Path(__file__).read_bytes()),"summary_sha256":sha((out/"summary.json").read_bytes()),"figures":["density-sensitivity.pdf","host-comparison.pdf"]})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--prepare",action="store_true");p.add_argument("--rebootstrap",action="store_true");p.add_argument("--replicates",type=int,default=10000);p.add_argument("--plot",action="store_true");p.add_argument("--lock",type=Path)
    a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    if a.prepare:prepare(a.root.resolve(),out)
    products,panel,subset,clusters,provenance=verified_inputs(out)
    fingerprint=sha(json.dumps({"schema":BOOTSTRAP_SCHEMA,"replicates":a.replicates,"seed":20260716,"products":provenance["checks"],"clusters":provenance["cluster_provenance"]["clusters_tsv_sha256"]},sort_keys=True).encode())
    bootstrap_path=out/"bootstrap-draws.csv.gz";bootstrap_meta=out/"bootstrap-provenance.json"
    if not bootstrap_path.exists() or a.rebootstrap:
        lock=(a.lock or out/"analysis.lock").resolve();lock.parent.mkdir(parents=True,exist_ok=True)
        with lock.open("a") as guard:
            print("Waiting for compute lock for bootstrap analysis",flush=True);fcntl.flock(guard,fcntl.LOCK_EX)
            draws=bootstrap(products,panel,clusters,a.replicates,20260716)
        save_csv(bootstrap_path,draws);save_json(bootstrap_meta,{"fingerprint":fingerprint,"draws_sha256":sha(bootstrap_path.read_bytes()),"script_sha256":sha(Path(__file__).read_bytes()),"schema":BOOTSTRAP_SCHEMA})
    else:
        meta=json.loads(bootstrap_meta.read_text())
        if meta["fingerprint"]!=fingerprint or meta["draws_sha256"]!=sha(bootstrap_path.read_bytes()):raise ValueError("Bootstrap cache does not match current frozen inputs/options")
        draws=read_csv(bootstrap_path)
    result=analyze(out,products,panel,subset,clusters,provenance,draws)
    if a.plot:plots(out)
    (out/"source").mkdir(exist_ok=True);shutil.copyfile(__file__,out/"source/publication_density.py")
    print(json.dumps({"rows":result["total_rows"],"clusters":result["clusters"],"human_fits":{k:v["fit"] for k,v in result["human"].items()}},indent=2),flush=True)


if __name__=="__main__":main()
