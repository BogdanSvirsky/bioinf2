#!/usr/bin/env python3
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VCF_PATH = ROOT / "student_35.vcf"
VEP_PATH = ROOT / "vep.json"
GNOMAD_PATH = ROOT / "gnomad.json"
CLINVAR_PATH = ROOT / "clinvar.json"
TABLE_PATH = ROOT / "annotated.tsv"
VEP_TABLE_PATH = ROOT / "vep_transcripts.tsv"
SETTINGS_PATH = ROOT / "annotation_settings.json"
USER_AGENT = "student-35-clinical-annotation/1.0"

VEP_PARAMS = {
    "symbol": 1,
    "canonical": 1,
    "mane": 1,
    "hgvs": 1,
    "numbers": 1,
    "sift": 1,
    "polyphen": 1,
    "af_gnomade": 1,
    "af_gnomadg": 1,
    "clin_sig": 1,
    "clin_sig_allele": 1,
}

CONSEQUENCE_RANK = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "TFBS_fusion",
    "TFBS_translocation",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_fusion",
    "regulatory_region_variant",
    "regulatory_region_truncation",
    "intergenic_variant",
    "sequence_variant",
]
RANK = {term: i for i, term in enumerate(CONSEQUENCE_RANK)}


def read_vcf():
    with VCF_PATH.open() as handle:
        rows = list(csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    return [
        {"chrom": r[0], "pos": int(r[1]), "id": r[2], "ref": r[3], "alt": r[4]}
        for r in rows
    ]


def request(url, payload=None, attempts=6, timeout=120):
    body = None if payload is None else json.dumps(payload).encode()
    headers = {"Accept": "application/json", "User-Agent": USER_AGENT}
    if body is not None:
        headers["Content-Type"] = "application/json"
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, data=body, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
            if attempt == attempts - 1:
                raise
            delay = min(2 ** attempt * 3, 30)
            print(f"Retry {attempt + 1}/{attempts - 1} after {type(error).__name__}: {url}")
            time.sleep(delay)


def fetch_vep(variants):
    if VEP_PATH.exists():
        return json.loads(VEP_PATH.read_text())
    results = []
    for start in range(0, len(variants), 25):
        ids = [v["id"] for v in variants[start:start + 25]]
        url = "https://rest.ensembl.org/vep/human/id?" + urllib.parse.urlencode(VEP_PARAMS)
        chunk = request(url, {"ids": ids})
        if not isinstance(chunk, list):
            raise RuntimeError(f"Unexpected VEP response: {chunk}")
        results.extend(chunk)
        print(f"VEP: {len(results)}/{len(variants)}")
    VEP_PATH.write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n")
    return results


def fetch_gnomad(variants):
    if GNOMAD_PATH.exists():
        return json.loads(GNOMAD_PATH.read_text())
    all_data = {}
    raw_chunks = []
    fields = "variant_id rsid chrom pos ref alt exome { ac an af } genome { ac an af }"
    for start in range(0, len(variants), 10):
        chunk = variants[start:start + 10]
        declarations = ", ".join(f"$v{i}: String!" for i in range(len(chunk)))
        selections = "\n".join(
            f'v{i}: variant(variantId: $v{i}, dataset: $dataset) {{ {fields} }}'
            for i in range(len(chunk))
        )
        query = f"query Variants({declarations}, $dataset: DatasetId!) {{\n{selections}\n}}"
        variables = {f"v{i}": f'{v["chrom"]}-{v["pos"]}-{v["ref"]}-{v["alt"]}' for i, v in enumerate(chunk)}
        variables["dataset"] = "gnomad_r4"
        response = request("https://gnomad.broadinstitute.org/api", {"query": query, "variables": variables})
        raw_chunks.append(response)
        if response.get("errors") and not response.get("data"):
            raise RuntimeError(f"gnomAD errors: {response['errors']}")
        for i, variant in enumerate(chunk):
            key = variant_key(variant)
            all_data[key] = response["data"][f"v{i}"]
        print(f"gnomAD: {len(all_data)}/{len(variants)}")
    output = {"dataset": "gnomad_r4", "variants": all_data, "raw_responses": raw_chunks}
    GNOMAD_PATH.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    return output


def fetch_clinvar(variants, vep_results):
    if CLINVAR_PATH.exists():
        return json.loads(CLINVAR_PATH.read_text())
    accessions = sorted({
        accession
        for vep in vep_results
        for colocated in vep.get("colocated_variants", [])
        for accession in colocated.get("var_synonyms", {}).get("ClinVar", [])
        if accession.startswith("VCV")
    })
    term = " OR ".join([v["id"] for v in variants] + accessions)
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    search_url = base + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "clinvar", "term": term, "retmode": "json", "retmax": 500}
    )
    search = request(search_url)
    ids = search["esearchresult"]["idlist"]
    summary = {"result": {"uids": []}}
    for start in range(0, len(ids), 100):
        chunk = ids[start:start + 100]
        summary_url = base + "esummary.fcgi?" + urllib.parse.urlencode(
            {"db": "clinvar", "id": ",".join(chunk), "retmode": "json"}
        )
        part = request(summary_url)
        if start == 0:
            summary = part
        else:
            summary["result"]["uids"].extend(part["result"]["uids"])
            for uid in part["result"]["uids"]:
                summary["result"][uid] = part["result"][uid]
    output = {"esearch": search, "esummary": summary}
    CLINVAR_PATH.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    print(f"ClinVar: {len(ids)} records")
    return output


def variant_key(variant):
    return f'{variant["chrom"]}:{variant["pos"]}:{variant["ref"]}>{variant["alt"]}'


def transcript_sort_key(transcript):
    terms = transcript.get("consequence_terms", [])
    severity = min((RANK.get(term, len(RANK)) for term in terms), default=len(RANK))
    return (
        severity,
        0 if transcript.get("mane_select") else 1,
        0 if transcript.get("mane_plus_clinical") else 1,
        0 if transcript.get("canonical") == 1 else 1,
        0 if transcript.get("biotype") == "protein_coding" else 1,
        transcript.get("transcript_id", ""),
    )


def select_transcript(vep, alt):
    transcripts = [
        t for t in vep.get("transcript_consequences", [])
        if t.get("variant_allele") == alt
    ]
    if not transcripts:
        raise RuntimeError(f'No VEP transcript for {vep.get("id")} allele {alt}')
    return min(transcripts, key=transcript_sort_key)


def clinvar_records(clinvar, variant):
    records = []
    result = clinvar["esummary"]["result"]
    expected_suffix = f':{variant["pos"] - 1}:{variant["ref"]}:{variant["alt"]}'
    for uid in result.get("uids", []):
        record = result[uid]
        for measure in record.get("variation_set", []):
            rsids = {
                "rs" + xref["db_id"]
                for xref in measure.get("variation_xrefs", [])
                if xref.get("db_source") == "dbSNP"
            }
            if variant["id"] in rsids or measure.get("canonical_spdi", "").endswith(expected_suffix):
                records.append(record)
                break
    return records


def format_clinvar(records, variant):
    matching = []
    for record in records:
        for measure in record.get("variation_set", []):
            spdi = measure.get("canonical_spdi", "")
            expected_suffix = f':{variant["pos"] - 1}:{variant["ref"]}:{variant["alt"]}'
            if spdi.endswith(expected_suffix):
                matching.append(record)
                break
    if not matching:
        return "нет", "—", "—"
    interpretations = []
    accessions = []
    for record in matching:
        accessions.append(record.get("accession_version", record.get("accession", "")))
        parts = []
        for field in ("germline_classification", "clinical_impact_classification", "oncogenicity_classification"):
            classification = record.get(field, {})
            description = classification.get("description", "")
            if description:
                review = classification.get("review_status", "")
                parts.append(f"{description} ({review})" if review else description)
        interpretations.append("; ".join(parts) if parts else "интерпретация не указана")
    return "да", "; ".join(dict.fromkeys(interpretations)), "; ".join(accessions)


def format_af(group):
    if not group or group.get("af") is None:
        return "NA"
    return f'{group["af"]:.8g}'


def build_outputs(variants, vep_results, gnomad, clinvar):
    vep_by_id = {}
    for item in vep_results:
        vep_by_id.setdefault(item["id"], []).append(item)
    rows = []
    transcript_rows = []
    for variant in variants:
        key = variant_key(variant)
        matching_vep = [
            item for item in vep_by_id.get(variant["id"], [])
            if item.get("assembly_name") == "GRCh38"
            and item.get("seq_region_name") == variant["chrom"]
            and item.get("start") == variant["pos"]
        ]
        if len(matching_vep) != 1:
            raise RuntimeError(f"Expected one primary GRCh38 VEP record for {key}, found {len(matching_vep)}")
        vep = matching_vep[0]
        selected = select_transcript(vep, variant["alt"])
        terms = sorted(selected["consequence_terms"], key=lambda term: RANK.get(term, len(RANK)))
        records = clinvar_records(clinvar, variant)
        present, interpretation, accessions = format_clinvar(records, variant)
        g = gnomad["variants"].get(key)
        rows.append({
            "variant": key,
            "source_id": variant["id"],
            "gene": selected.get("gene_symbol", selected.get("gene_id", "NA")),
            "vep_consequence": ",".join(terms),
            "transcript": selected.get("transcript_id", "NA"),
            "mane_select": selected.get("mane_select", "NA"),
            "canonical": selected.get("canonical", "NA"),
            "hgvsc": selected.get("hgvsc", "NA"),
            "hgvsp": selected.get("hgvsp", "NA"),
            "sift": selected.get("sift_prediction", "NA"),
            "polyphen": selected.get("polyphen_prediction", "NA"),
            "rsid": variant["id"],
            "gnomad_exome_af": format_af((g or {}).get("exome")),
            "gnomad_genome_af": format_af((g or {}).get("genome")),
            "clinvar_present": present,
            "clinvar_interpretation": interpretation,
            "clinvar_accessions": accessions,
        })
        for transcript in vep.get("transcript_consequences", []):
            if transcript.get("variant_allele") != variant["alt"]:
                continue
            transcript_rows.append({
                "variant": key,
                "source_id": variant["id"],
                "gene": transcript.get("gene_symbol", transcript.get("gene_id", "NA")),
                "transcript": transcript.get("transcript_id", "NA"),
                "consequences": ",".join(transcript.get("consequence_terms", [])),
                "impact": transcript.get("impact", "NA"),
                "mane_select": transcript.get("mane_select", "NA"),
                "canonical": transcript.get("canonical", "NA"),
                "hgvsc": transcript.get("hgvsc", "NA"),
                "hgvsp": transcript.get("hgvsp", "NA"),
                "sift": transcript.get("sift_prediction", "NA"),
                "polyphen": transcript.get("polyphen_prediction", "NA"),
            })
    with TABLE_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with VEP_TABLE_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(transcript_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(transcript_rows)
    return rows


def write_settings():
    settings = {
        "date": date.today().isoformat(),
        "input": VCF_PATH.name,
        "assembly": "GRCh38",
        "vep": {
            "resource": "Ensembl VEP REST API",
            "endpoint": "https://rest.ensembl.org/vep/human/id",
            "method": "POST, rsID identifiers",
            "parameters": VEP_PARAMS,
            "raw_output": VEP_PATH.name,
            "all_transcripts_output": VEP_TABLE_PATH.name,
        },
        "clinvar": {
            "resource": "NCBI ClinVar E-utilities",
            "search": "db=clinvar; input rsIDs and VCV accessions from VEP combined with OR; retmax=500; retmode=json",
            "summary": "db=clinvar; ESummary for every search hit; retmode=json",
            "allele_matching": "GRCh38 canonical SPDI must match input position, REF and ALT",
            "raw_output": CLINVAR_PATH.name,
        },
        "gnomad": {
            "resource": "gnomAD GraphQL API",
            "endpoint": "https://gnomad.broadinstitute.org/api",
            "dataset": "gnomad_r4",
            "variant_id_format": "CHROM-POS-REF-ALT (GRCh38)",
            "fields": "variant_id, rsid, chrom, pos, ref, alt, exome/genome AC, AN, AF",
            "raw_output": GNOMAD_PATH.name,
        },
        "transcript_selection_rule": [
            "Only transcript consequences for the exact input ALT allele",
            "Most severe consequence according to the Ensembl consequence ranking",
            "For equal severity prefer MANE Select, then MANE Plus Clinical, then canonical, then protein_coding, then lexicographically smallest transcript ID",
        ],
        "summary_output": TABLE_PATH.name,
    }
    SETTINGS_PATH.write_text(json.dumps(settings, ensure_ascii=False, indent=2) + "\n")


def main():
    variants = read_vcf()
    if len(variants) != 50:
        raise RuntimeError(f"Expected 50 variants, found {len(variants)}")
    vep = fetch_vep(variants)
    gnomad = fetch_gnomad(variants)
    clinvar = fetch_clinvar(variants, vep)
    rows = build_outputs(variants, vep, gnomad, clinvar)
    write_settings()
    print(f"Wrote {len(rows)} rows to {TABLE_PATH.name}")


if __name__ == "__main__":
    main()
