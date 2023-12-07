import os
import pathlib
import subprocess

this_file_path = os.path.dirname(os.path.realpath(__file__))
print(this_file_path)

should_build = False
run_no_gui = True

# if should_build:
#     subprocess.run("cmake -B build", cwd=this_file_path)
#     subprocess.run("cmake --build build", cwd=this_file_path)
if run_no_gui:
    subprocess.run("./bin/SPHSimulator --no-gui", cwd=this_file_path)