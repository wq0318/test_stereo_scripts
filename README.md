# test_stereo_scripts

This directory is the script root copied into the docker image as `/test_stereo_scripts`.

## Layout

- `analysis/`
  - `downstream_analysis.sh`: stable entrypoint used by WDL
  - `downstream_analysis_impl.sh`: current downstream implementation
  - `atac_saturation.pl`: saturation and cutoff helper
  - `draw_atac_qc.R`: QC plotting and spatial filtering
  - `signac_analysis.R`: Signac downstream analysis
- `reporting/`
  - `parse_barcode_mapping_log.py`: normalize ST_BarcodeMap log into metrics CSV
  - `parse_mapping_summary.py`: normalize chromap summary CSV into metrics CSV
  - `parse_sam2frag_report.py`: normalize PISA sam2bam report into metrics CSV
  - `parse_merge_report.py`: normalize merged fragment stats into metrics CSV
  - `parse_downstream_analysis_outputs.py`: normalize downstream qc/signac outputs into metrics CSV
  - `render_final_report.py`: generate the final HTML report for each bin size
- `tests/`
  - `test_barcode_mapping.sh`: standalone barcodeMapping test
  - `test_mapping.sh`: standalone chromap mapping test
  - `test_sam2frag.sh`: standalone sam2bam+bam2frag test
  - `test_merge_frag.sh`: standalone merge fragment test
  - `test_downstream_analysis.sh`: standalone downstream analysis test
  - `test_render_final_report.sh`: final html rendering test
  - `stepwise_test_template.sh`: fill in paths and run step by step on the cloud platform

## Naming

- raw task outputs keep the original sample-oriented names
- normalized reports use `*.metrics.csv`
- downstream entry scripts use verb-oriented names such as `downstream_analysis.sh`

## Output Layout

- `01.barcode_mapping/`
- `02.mapping/`
- `03.sam2frag/`
- `04.merge/`
- `05.bin_square/`
- `06.analysis/bin50` and `06.analysis/bin100`
- `07.report/bin50` and `07.report/bin100`
