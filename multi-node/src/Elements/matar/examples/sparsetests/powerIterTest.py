import subprocess, os

def powerIterations():
    subprocess.call('rm ../results/powerIterations.csv', shell=True)
    Ns = [1e3, 1e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7]
    subprocess.call('echo N, Dense Time(s), Sparse Time(s) >> ../results/powerIterations.csv', shell=True)
    for n in Ns:
        for i in range(3):
            print("n: ", n, "loop: ", i, " of ", 3)
            subprocess.call('./examples/sparsetests/powerIter ' + n  + ' >> results/powerIterations.csv', shell=True);



powerIterations()
