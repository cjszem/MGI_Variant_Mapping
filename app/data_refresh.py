import yaml
from data_util import download_file, parse_mondo_obo, parse_ensembl_gff3


# Load Config
with open('app/config.yaml') as f:
    config = yaml.safe_load(f)


def refresh_mondo():
    url = config['data']['mondo_obo']
    output_path = config['paths']['mondo_term_map']

    download_file(url, output_path)
    parse_mondo_obo(output_path)


def refresh_mgi():
    download_file(config['data']['mgi_allele'])
    download_file(config['data']['mgi_homology'])
    download_file(config['data']['mgi_disease'])
    download_file(config['data']['mgi_phenotype'])


def refresh_clinvar():
    url_vcf = config['data']['clinvar_vcf']
    url_tbi = config['data']['clinvar_tbi']
    out_vcf = config['paths']['clinvar_vcf']
    out_tbi = config['paths']['clinvar_tbi']

    download_file(url_vcf, out_vcf)
    download_file(url_tbi, out_tbi)


def refresh_ensembl():
    url_hum = config['data']['ensembl_hum_gff3']
    url_mus = config['data']['ensembl_mus_gff3']
    out_hum = config['paths']['ensembl_hum_gff3']
    out_mus = config['paths']['ensembl_mus_gff3']

    download_file(url_hum, out_hum)
    download_file(url_mus, out_mus)


def refresh_data():
    refresh_mondo()
    refresh_mgi()
    refresh_clinvar()
    refresh_ensembl()