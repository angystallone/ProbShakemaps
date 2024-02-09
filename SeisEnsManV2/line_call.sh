mainFolder=$(pwd)
python write_json_file.py
python run_ens.py --cfg $mainFolder/input/main.config --event $mainFolder/input/event_stat.json --nb_scen 1000
