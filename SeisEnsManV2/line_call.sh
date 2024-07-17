mainFolder=$(pwd)

python write_json_file.py

exit_status=$?
if [ $exit_status -ne 0 ]; then
    echo "Error: write_json_file.py exited with non-zero status code: $exit_status"
    exit $exit_status
fi

python run_ens.py --cfg $mainFolder/input/main.config --event $mainFolder/input/event_stat.json --nb_scen 1000 --angles 318 89 -179

# norcia 151 47 -89
# turchia 318 89 -179

