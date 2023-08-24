import subprocess, os

device = "amd"

def runParN():
    subprocess.call('rm ../results/results_par_' + device + '.csv', shell=True)

    Ns = [500,1000, 2000, 4000, 6000, 8000, 10000]

    Ns = [str(n) for n in Ns] 

    subprocess.call('echo N, Prob, K, t1, t2, t3, total_time, distance >> ../results/results_par_' + device + '.csv', shell=True)
    for n in Ns:
        for i in range(3):
            print("n:", n , "loop:", i, " of 3")
            subprocess.call('./examples/watt-graph/test_kokkos_floyd ' + n + ' 0 6 >> ../results/results_par_' + device + '.csv', shell=True)


def runDiffP():
    subprocess.call('rm ../results/results_dif_p.csv', shell=True)

    Ps = [0.0001, 0.005, 0.01, 0.02, 0.05, 0.1]

    Ps = [str(p) for p in Ps] 

    subprocess.call('echo N, Prob, K, t1, t2, t3, total_time, distance >> ../results/results_dif_p.csv', shell=True)
    for p in Ps:
        for i in range(3):
            print("p:", p , "loop:", i, " of 3")
            subprocess.call('./examples/watt-graph/test_kokkos_floyd 5000 '  + p + ' 6 >> ../results/results_dif_p.csv', shell=True)



def runSerN():
    subprocess.call('rm  ../results/results_ser_' + device + '.csv', shell=True)

    Ns = [500,1000, 2000, 2500,  3000, 3500, 4000]

    Ns = [str(n) for n in Ns] 

    subprocess.call('echo N, Prob, K, t1, t2, t3, total_time, distance >> ../results/results_ser_' + device + '.csv', shell=True)
    for n in Ns:
        for i in range(3):
            print("n:", n , "loop:", i, " of 3")
            subprocess.call('./examples/watt-graph/test_floyd ' + n + ' 0 6 >>  ../results/results_ser_' + device + '.csv' , shell=True)



def runCpuN():

    subprocess.call('rm ../results/results_cpu.csv', shell=True)

    Ns = [500,1000, 2000, 2500,  3000, 3500, 4000, 6000]

    Ns = [str(n) for n in Ns] 

    subprocess.call('echo N, Prob, K, t1, t2, t3, total_time, distance >> ../results/results_cpu.csv', shell=True)
    for n in Ns:
        for i in range(3):
            print("n:", n , "loop:", i, " of 3")
            subprocess.call('./examples/watt-graph/test_kokkos_floyd ' + n + ' 0 6 >> ../results/results_cpu.csv', shell=True)


runCpu()
