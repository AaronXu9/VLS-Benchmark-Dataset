from chembl_webresource_client.new_client import new_client
import pandas as pd
from time import time
activity = new_client.activity
assay = new_client.assay

start_time = time()
acts = activity.filter(
    target_organism='Homo sapiens',
    standard_type__in=['IC50', 'Ki', 'Kd', 'EC50'],
    target_components__accession="P0053"
)
acts = acts[:20]
assay_ids = {a['assay_chembl_id'] for a in acts if a.get('assay_chembl_id')}
assays = assay.filter(assay_chembl_id__in=list(assay_ids)).only(
    'assay_chembl_id', 'assay_organism', 'assay_type', 'assay_tissue', 'assay_cell_type', 'assay_category', 'assay_strain', 'assay_cell_type', 'assay_tax_id', 'assay_tissue', 'assay_subcellular_fraction'
)
end_time = time()
print(end_time-start_time)

assays_df = pd.DataFrame(assays)
# then merge on 'assay_chembl_id'
print(assays_df)