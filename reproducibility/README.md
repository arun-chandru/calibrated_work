# numerical reproduction

These scripts recompute the finite numerical certificates supporting the manuscript. They do not compile a manuscript or generate a PDF. The analytic derivations are in the manuscript; executable output does not replace those derivations.

All five Python files are unchanged copies of the verified source scripts. They require only Python's standard library. This package was tested with Python 3.12.14. Keep the scripts and three supplied reference JSON files together in this directory.

Run the commands sequentially from this directory. Keep Python assertions enabled: do not use `-O` or set `PYTHONOPTIMIZE` to a positive value. `-B` prevents creation of bytecode cache files.

```sh
python -c "from pathlib import Path; Path('outputs').mkdir(exist_ok=True)"
python -B verify_uv_smooth_calibrated_intervals.py --depth 12 --output outputs/smooth_reproduced.json
python -B verify_uv_preparation_comparison.py --depth 12 --output outputs/preparations_reproduced.json
python -B verify_uv_bounded_design_intervals.py --depth 12 --output outputs/design_reproduced.json
python -B audit_uv_bounded_design.py
```

Depth 12 means 4096 integration cells. The bounded-design command recomputes all eight candidates. The last command checks the supplied smooth/design records and freshly rebuilds the selected candidate; it creates `PHY_UV_BOUNDED_DESIGN_INDEPENDENT_QA_2026_10_06.json` beside the scripts. The smooth, preparation, and design generation commands compute their results without reading the supplied reference JSON files.

The three reference files for comparison are:

- `PHY_UV_SMOOTH_CALIBRATED_CERTIFICATE_2026_10_06.json`
- `PHY_UV_PREPARATION_COMPARISON_CERTIFICATE_2026_10_06.json`
- `PHY_UV_BOUNDED_DESIGN_CERTIFICATE_2026_10_06.json`

Compare the complete parsed results with this optional standard-library command. It exits unsuccessfully on any difference:

```sh
python -B -c "import json; from pathlib import Path; pairs=[('PHY_UV_SMOOTH_CALIBRATED_CERTIFICATE_2026_10_06.json','outputs/smooth_reproduced.json'),('PHY_UV_PREPARATION_COMPARISON_CERTIFICATE_2026_10_06.json','outputs/preparations_reproduced.json'),('PHY_UV_BOUNDED_DESIGN_CERTIFICATE_2026_10_06.json','outputs/design_reproduced.json')]; assert all(json.loads(Path(a).read_text(encoding='utf-8')) == json.loads(Path(b).read_text(encoding='utf-8')) for a,b in pairs), 'Reproduction differs from a reference certificate'; print('All three complete parsed certificates match')"
```

Parsed equality is the portable comparison: operating-system newline conventions can change byte hashes. The preparation output covers six declared preparations. The design selection concerns the best certified detection guarantee within the declared eight-element family; it does not establish a global optimum or actual-probability superiority.
