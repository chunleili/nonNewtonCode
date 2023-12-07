from subprocess import call
# cmake --preset vs2022
# cmake --build --preset vs2022-Rel
# ./bin/SPHSimulator --scene-file ../data/MyScenes/ramp1.json --no-initial-pause
# ./bin/SPHSimulator --scene-file ../data/MyScenes/strawberry_honey.json --no-initial-pause
call(["cmake", "--preset", "vs2022"])
call(["cmake", "--build", "--preset", "vs2022-Rel"])
call(["./bin/SPHSimulator", "--scene-file", "../data/MyScenes/strawberry_honey.json", "--no-initial-pause"])