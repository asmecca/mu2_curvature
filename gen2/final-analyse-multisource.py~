#!/usr/bin/python3
# coding=utf-8
"""
Analyse the contents of the mu^2 computation
"""

import os
import re
import sys

import numpy as np

N_gamma = 16
N_time = int(sys.argv[2])

a_inv_gev = 6.079 #in GeV
T = float(sys.argv[4]) #in GeV
T_lat = float(float(T)/float(a_inv_gev))
mu = float(sys.argv[3])#0.281 #in GeV
mu_lat = float(float(mu)/float(a_inv_gev))


def extract_run_info(dirname):
    """Extract the run information from the files in the folder

    Keyword arguments:
    dirname -- directory data is stored

    Output:
     - name of the simulation
     - list of mass id's
     - list of configuration id's
    """

    root_name = ""
    configs = []
    masses = []
    sources = []

    for file in os.listdir(dirname):
        m = re.match(r'(.*?)n(\d+)\.(one-all\.s(\d+)|all-all)\.m(\d+)', file)

        if not m:
            continue

        if root_name == "":
            root_name = m.group(1)
        elif root_name != m.group(1):
            raise RuntimeError("Inconsistent run names")

        configs.append(int(m.group(2)))

        if m.group(4):
            sources.append(int(m.group(4)))

        masses.append(int(m.group(5)))


    return (root_name, sorted(set(masses)), np.array(sorted(set(configs))), sorted(set(sources)))


def one_all_data(dirname, run_name, masses, cfgs, srcs):
    """Read the one-all contribution data

    Keyword arguemnts:
    dirname -- directory data is stored
    run_name -- name of the simulation
    masses -- list of mass ID's to read
    cfgs -- list of configuration numbers to read
    srcs -- list of one-all sources

    Output:
    One-All data in a matrix with dimensions:
    Nmasses x Nconfigs x Ngamma x Ntimeslices x Nterms

    Averaged over sources
    """

    data = np.zeros(
        (len(masses), len(cfgs), N_gamma, N_time, 5), dtype=complex)

    for im, m in enumerate(masses):
        for icfg, cfg in enumerate(cfgs):
            for src in srcs:
                tmp_dat = np.loadtxt("{}/{}n{}.one-all.s{}.m{}".format(
                    dirname, run_name, cfg, src, m)).reshape(N_time, N_gamma, 12)

                tmp_dat = tmp_dat[:, :, 2::2] + 1j * tmp_dat[:, :, 3::2]
                data[im, icfg] += np.transpose(tmp_dat, (1, 0, 2))

            data[im, icfg] /= len(srcs)

    return data


def all_all_statistics(dirname, run_name, masses, cfgs):
    """Read the all-all contribution data and take the statistical average of
    the noise vectors

    Keyword arguemnts:
    dirname -- directory data is stored
    run_name -- name of the simulation
    masses -- list of mass ID's to read
    cfgs -- list of configuration numbers to read

    Output:
    All-All statistical average
    Nmasses x Nconfigs x Nterms x {avg, std}
    """

    stats = np.zeros((len(masses), len(cfgs), 4), dtype=complex)

    for im, m in enumerate(masses):
        for icfg, cfg in enumerate(cfgs):

            # Compute the 3 traces estimated with noise vectors
            dat = np.loadtxt(
                "{}/{}n{}.all-all.m{}".format(dirname, run_name, cfg, m))

            n = len(dat)

            dat = dat[:, ::2] + 1j * dat[:, 1::2]
            stats[im, icfg, 0:3] = np.mean(dat, axis=0)

            # Compute the tr[dot{D} D^{-1}]^2 term
            dat = np.power(dat[:, 0], 2.0)

            tmp_sqr = 1.0 / (n * (
                n - 1.0)) * ((n * stats[im, icfg, 0])**2 - np.sum(dat))

            stats[im, icfg, 3] = tmp_sqr

    return stats


def compute_contributions(oa_data, aa_data, bootstrap_samples=10000):
    """Compute the different terms of the mu^2 calculations

    They are:
    0. Re < oa[4] > (the raw contribution)
    1. +4.0 Re < oa[0] >
    2. -2.0 Re < oa[1] >
    3. -2.0 Re < oa[2] >
    4. -8.0 < Im(oa[3]) Im(aa[0]) >
    5. +2.0 < oa[4] * (aa[1] - aa[2] + 2*aa[3]) - <oa[4]> <aa[1] - aa[2] + 2*aa[3]> >
    6. Sum of the above

    Configuration averages and errors are estimated with bootstrap

    Keyword arguemnts:
    oa_data -- the one-all data
    aa_data -- the all-all data
    bootstrap_samples -- number of bootstrap samples to use (default=10000)

    Output:
    The different terms in a matrix of dimensions
    N_masses x N_gamma x N_time x N_terms(6) x {avg, std}
    """

    num_masses = aa_data.shape[0]
    N = aa_data.shape[1]

    res = np.zeros((num_masses, N_gamma, N_time, 8, 2))
    vec = np.zeros((num_masses, 1, N_time, 8, 2))
    axial = np.zeros((num_masses, 1, N_time, 8, 2))
    g_ii = np.zeros((num_masses, N_gamma, N_time, 1, 2))
    g_ii_vec = np.zeros((num_masses, 1, N_time, 1, 2))
    g_ii_axial = np.zeros((num_masses, 1, N_time, 1, 2))
    disc = np.zeros((num_masses, N_gamma, N_time, 1, 2))
    disc_vec = np.zeros((num_masses, 1, N_time, 1, 2))
    disc_axial = np.zeros((num_masses, 1, N_time, 1, 2))    
    
    for im in range(num_masses):
        print("THIS IS THE NUMB OF SAMPLES: ",bootstrap_samples)
        samples = [np.random.choice(N, N) for s in range(bootstrap_samples)]
        vals = np.zeros((bootstrap_samples, N_gamma, N_time, 8))
        vec_vals = np.zeros((bootstrap_samples, 1, N_time, 8))
        axial_vals = np.zeros((bootstrap_samples, 1, N_time, 8))
        g_ii_vals = np.zeros((bootstrap_samples, N_gamma, N_time, 1))
        g_ii_vals_vec = np.zeros((bootstrap_samples, 1, N_time, 1))
        g_ii_vals_axial = np.zeros((bootstrap_samples, 1, N_time, 1))
        disc_vals = np.zeros((bootstrap_samples, N_gamma, N_time, 1))
        disc_vals_vec = np.zeros((bootstrap_samples, 1, N_time, 1))
        disc_vals_axial = np.zeros((bootstrap_samples, 1, N_time, 1))

        for s_id, s in enumerate(samples):

            aa_sample = aa_data[im, s]
            oa_sample = oa_data[im, s]

            # Zeroth term
            vals[s_id, :, :, 0] = np.mean(
                np.real(oa_sample[:, :, :, 4]), axis=0)

            # First term
            vals[s_id, :, :, 1] = 4 * np.mean(
                np.real(oa_sample[:, :, :, 0]), axis=0)

            # Second term
            vals[s_id, :, :, 2] = -2 * np.mean(
                np.real(oa_sample[:, :, :, 1]), axis=0)
            
            # Third term
            vals[s_id, :, :, 3] = -2 * np.mean(
                np.real(oa_sample[:, :, :, 2]), axis=0)

            # Fourth term
            vals[s_id, :, :, 4] = -8 * np.mean(
                np.multiply(
                    np.imag(np.transpose(oa_sample[:, :, :, 3], (1, 2, 0))),
                    np.imag(aa_sample[:, 0])),
                axis=2)

            # Fifth term
            oa_sample_avg = np.mean(oa_sample[:, :, :, 4], axis=0)
            aa_sample_avg = np.mean(aa_sample[:, 1] - aa_sample[:, 2] +
                                    2 * aa_sample[:, 3])

            vals[s_id, :, :, 5] = 2 * np.real(
                np.mean(
                    np.multiply(
                        np.transpose(oa_sample[:, :, :, 4] - oa_sample_avg, (
                            1, 2, 0)),
                        np.real(aa_sample[:, 1] - aa_sample[:, 2] +
                                2 * aa_sample[:, 3]) - aa_sample_avg),
                    axis=2))

            # Sixth term (sum)
            vals[s_id, :, :, 6] = np.sum(vals[s_id, :, :, 1:6], axis=2)

            # Seventh term ( O(1)+O(mu^2) i.e. sum of G and G'' )
            vals[s_id, :, :, 7] = np.add(vals[s_id, :, :, 6]*0.5*((mu_lat)**2), vals[s_id, :, :, 0])

            #computing vec
            vec_vals[s_id, 0, :, :] = np.add(np.add(vals[s_id, 2, :, :],vals[s_id, 3, :, :]),vals[s_id, 4 , :, :])/3.0

            #computing axial
            axial_vals[s_id, 0, :, :] = np.add(np.add(vals[s_id, 7, :, :],vals[s_id, 8, :, :]),vals[s_id, 9 , :, :])/3.0

            #computing connected
            g_ii_vals[s_id, :, :, 0] = np.add(np.add(vals[s_id, :, :, 1],vals[s_id, :, :, 2]),vals[s_id, :, :, 3])

            #computing connected vector
            g_ii_vals_vec[s_id, 0, :, 0] = np.add( np.add(
                np.add( np.add( vals[s_id, 2, :, 1], vals[s_id, 3, :, 1]), vals[s_id, 4, :, 1] )/3.0 ,
                np.add( np.add( vals[s_id, 2, :, 2], vals[s_id, 3, :, 2]), vals[s_id, 4, :, 2] )/3.0 ),
                np.add( np.add( vals[s_id, 2, :, 3], vals[s_id, 3, :, 3]), vals[s_id, 4, :, 3] )/3.0)

            #computing connected axial
            g_ii_vals_axial[s_id, 0, :, 0] = np.add( np.add(
                np.add( np.add( vals[s_id, 7, :, 1], vals[s_id, 8, :, 1]), vals[s_id, 9, :, 1] )/3.0 ,
                np.add( np.add( vals[s_id, 7, :, 2], vals[s_id, 8, :, 2]), vals[s_id, 9, :, 2] )/3.0 ),
                np.add( np.add( vals[s_id, 7, :, 3], vals[s_id, 8, :, 3]), vals[s_id, 9, :, 3] )/3.0)
                        
            #computing disconnected
            disc_vals[s_id, :, :, 0] = np.add(vals[s_id, :, :, 4],vals[s_id, :, :, 5])

            #computing connected vector            
            disc_vals_vec[s_id, 0, :, 0] = np.add(
                np.add( np.add( vals[s_id, 2, :, 4], vals[s_id, 3, :, 4] ), vals[s_id, 4, :, 4])/3.0,
                np.add( np.add( vals[s_id, 2, :, 5], vals[s_id, 3, :, 5] ), vals[s_id, 4, :, 5])/3.0 )

            #computing connected axial
            disc_vals_axial[s_id, 0, :, 0] = np.add(
                np.add( np.add( vals[s_id, 7, :, 4], vals[s_id, 8, :, 4] ), vals[s_id, 9, :, 4])/3.0,
                np.add( np.add( vals[s_id, 7, :, 5], vals[s_id, 8, :, 5] ), vals[s_id, 9, :, 5])/3.0 )

            #Saving each bootstrap sample in directory boot for vector and axial-vector correlators
            res_full_dirname = "{}/analysis_mu_{}/boot".format(sys.argv[1], mu)            
            if not os.path.exists(res_full_dirname):
                os.makedirs(res_full_dirname)
            for t in range(0, N_time):
                with open("{}/res.vector.dat".format(res_full_dirname),"a") as f:                    
                    np.savetxt( f ,
                               vec_vals[s_id, 0, t].reshape(1,8))                    
                f.close()
                with open("{}/res.axial.dat".format(res_full_dirname),"a") as f:
                    np.savetxt( f ,
                                axial_vals[s_id, 0, t].reshape(1,8))
                f.close()        
        
        # Bootstrap estimate for mean and error
        res[im, :, :, :, 0] = np.mean(vals, axis=0)
        res[im, :, :, :, 1] = np.std(vals, axis=0)
        vec[im, :, :, :, 0] = np.mean(vec_vals, axis=0)
        vec[im, :, :, :, 1] = np.std(vec_vals, axis=0)
        axial[im, :, :, :, 0] = np.mean(axial_vals, axis=0)
        axial[im, :, :, :, 1] = np.std(axial_vals, axis=0)
        g_ii[im, :, :, :, 0] = np.mean(g_ii_vals, axis=0)
        g_ii[im, :, :, :, 1] = np.std(g_ii_vals, axis=0)
        g_ii_vec[im, :, :, :, 0] = np.mean(g_ii_vals_vec, axis=0)
        g_ii_vec[im, :, :, :, 1] = np.std(g_ii_vals_vec, axis=0)
        g_ii_axial[im, :, :, :, 0] = np.mean(g_ii_vals_axial, axis=0)
        g_ii_axial[im, :, :, :, 1] = np.std(g_ii_vals_axial, axis=0)
        disc[im, :, :, :, 0] = np.mean(disc_vals, axis=0)
        disc[im, :, :, :, 1] = np.std(disc_vals, axis=0)
        disc_vec[im, :, :, :, 0] = np.mean(disc_vals_vec, axis=0)
        disc_vec[im, :, :, :, 1] = np.std(disc_vals_vec, axis=0)
        disc_axial[im, :, :, :, 0] = np.mean(disc_vals_axial, axis=0)
        disc_axial[im, :, :, :, 1] = np.std(disc_vals_axial, axis=0)

    return res,vec,axial,g_ii,g_ii_vec,g_ii_axial,disc,disc_vec,disc_axial


def main(args):
    """
    The main data analysis

    Keyword arguments:
    args -- List of command line arguments, assumed to be the directory in
            which the data output of the mu-squared code is stored

    Output:
    Stores the results in a folder args[1]/analysis
    """

    if not os.path.isdir(args[1]):
        raise RuntimeError(args[1] +
                           " is does not exist or is not a directory")

    num_samples = 2000

    if (len(args) > 5):
        num_samples = int(args[6])
        

    run_name, masses, cfgs, sources = extract_run_info(args[1])
    print(run_name, masses, cfgs, sources)

    aa_stats = all_all_statistics(args[1], run_name, masses, cfgs)
    oa_data = one_all_data(args[1], run_name, masses, cfgs, sources)

    res,vec,axial,g_ii,g_ii_vec,g_ii_axial,disc,disc_vec,disc_axial = compute_contributions(oa_data, aa_stats, num_samples)

    res_full_dirname = "{}/analysis_mu_{}/full".format(args[1], mu)
    if not os.path.exists(res_full_dirname):
        os.makedirs(res_full_dirname)

    for im in range(len(masses)):
        for ig in range(N_gamma):
            np.savetxt("{}/res.g{}.m{}.dat".format(res_full_dirname, ig, im),
                       res[im, ig, :, 6].reshape(N_time, 2))
        np.savetxt("{}/vec.m{}.dat".format(res_full_dirname, im),
                   vec[im, 0, :, 6].reshape(N_time, 2))
        np.savetxt("{}/axial.m{}.dat".format(res_full_dirname, im),
                   axial[im, 0, :, 6].reshape(N_time, 2))

    res_terms_dirname = "{}/analysis_mu_{}/terms".format(args[1], mu)
    if not os.path.exists(res_terms_dirname):
        os.makedirs(res_terms_dirname)

    for im in range(len(masses)):
        for ig in range(N_gamma):
            np.savetxt("{}/res.g{}.m{}.dat".format(res_terms_dirname, ig, im),
                       res[im, ig].reshape(N_time, 16))
        np.savetxt("{}/res_g_ii.m{}.dat".format(res_terms_dirname, im),
                       g_ii[im, 0, :, 0].reshape(N_time, 2))
        np.savetxt("{}/res_disc.m{}.dat".format(res_terms_dirname, im),
                   disc[im, 0, :, 0].reshape(N_time, 2))
        np.savetxt("{}/res.vector.m{}.dat".format(res_terms_dirname, im),
                   vec[im, 0].reshape(N_time, 16))
        np.savetxt("{}/res.axial.m{}.dat".format(res_terms_dirname, im),
                   axial[im, 0].reshape(N_time, 16))
        np.savetxt("{}/res.g_ii_vector.m{}.dat".format(res_terms_dirname, im),
                   g_ii_vec[im, 0, :, 0].reshape(N_time, 2))
        np.savetxt("{}/res.g_ii_axial.m{}.dat".format(res_terms_dirname, im),
                   g_ii_axial[im, 0, :, 0].reshape(N_time, 2))
        np.savetxt("{}/res.disc_vector.m{}.dat".format(res_terms_dirname, im),
                   g_ii_vec[im, 0, :, 0].reshape(N_time, 2))
        np.savetxt("{}/res.disc_axial.m{}.dat".format(res_terms_dirname, im),
                   g_ii_axial[im, 0, :, 0].reshape(N_time, 2))


if __name__ == "__main__":
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        print("\nProgram interrupted by user")
