import subprocess, os
import time 

device = "turing"

def matVec():
    f = '../results/sparse/matVec_' + device + '.csv'
    subprocess.call('rm ' + f, shell=True)

    Ns = [500,1000, 2000, 4000, 6000, 8000, 10000, 20000, 30000, 40000, 50000, 60000]

    Ns = [str(n) for n in Ns] 

    subprocess.call('echo N, Prob, K, t1, t2, t3, total_time, distance >>' + f , shell=True)
    for n in Ns:
        for i in range(3):
            print("n:", n , "loop:", i, " of 3")
            subprocess.call('./examples/sparsetests/matVec ' + n + ' >> ' + f, shell=True)
            time.sleep(0.5)

def spatVec():
    f = '../results/sparse/spatVec_' + device + '.csv'
    subprocess.call('rm ' + f, shell=True)

    Ns = [i*10000 for i in range(1,200)]

    Ns = [str(n) for n in Ns] 

    subprocess.call('echo N, Prob, K, t1, t2, t3, total_time, distance >>' + f , shell=True)
    for n in Ns:
        for i in range(3):
            print("n:", n , "loop:", i, " of 3")
            subprocess.call('./examples/sparsetests/spatVec ' + n + ' >> ' + f, shell=True)
            time.sleep(0.5)

def powerIter():
    f = '../results/sparse/powerIter_' + device + '.csv'
    subprocess.call('rm ' + f, shell=True)

    Ns = [100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000,12000,1300,14000,15000,16000,17000,18000,19000,20000]

    Ns = [str(n) for n in Ns] 

    subprocess.call('echo N, Prob, K, t1, t2, t3, total_time, distance >>' + f , shell=True)
    for n in Ns:
        for i in range(3):
            print("n:", n , "loop:", i, " of 3")
            subprocess.call('./examples/sparsetests/powerIter ' + n + ' >> ' + f, shell=True)
            time.sleep(0.5)



powerIter()
