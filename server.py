from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.encoders import jsonable_encoder

from app.main_util import (batch_hvar, batch_mvar, batch_score)
from app.processing_util import process_batch_query

app = FastAPI()

app.mount("/static", StaticFiles(directory="static"), name="static")

@app.get("/")
def serve_dashboard():
    return FileResponse("templates/dashboard.html")

class VariantInput(BaseModel):
    input: str

@app.post("/run_variants")
def run_variants(data: VariantInput):

    print(data.input)

    variants = process_batch_query(data.input)

    hum_gene_df, hum_prt_df, input_gene_df = batch_hvar(variants)
    mouse_gene_df, mouse_prt_df, gene_input_df = batch_mvar(input_gene_df)
    score_df = batch_score(hum_prt_df, mouse_prt_df, gene_input_df)

    # results = {"human_genes": hum_gene_df.to_dict(orient="records"),
    #            "human_proteins": hum_prt_df.to_dict(orient="records"),
    #            "mouse_genes": mouse_gene_df.to_dict(orient="records"),
    #            "mouse_proteins": mouse_prt_df.to_dict(orient="records"),
    #            "scores": score_df.to_dict(orient="records")}

    # safe_result = jsonable_encoder(results)

    html_results = {
        "human_genes": hum_gene_df.to_html(index=False),
        "human_proteins": hum_prt_df.to_html(index=False),
        "mouse_genes": mouse_gene_df.to_html(index=False),
        "mouse_proteins": mouse_prt_df.to_html(index=False),
        "scores": score_df.to_html(index=False)
    }

    return JSONResponse(content=html_results)
