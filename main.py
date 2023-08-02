import csv, json, math, requests, time

from load_scaffold_data import read_herb_data
from prep_data import prep_data
from prep_data import prep_data
from prep_data import prep_data
from prep_data import prep_data
from prep_data import prep_data
from prep_data import prep_data

import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

#DA_setting
API_Key = <API_Key>
API_Access_URL = <API_Access_URL>
headers = {'content-type':'application/json', 'X-Api-Key':API_Key}
proxies = {}
A = 1000  # weight for [one sidechain for one site] penalty terms
B = 1000  # weight for [weight equaling] penalty terms
C = 1  # weight for [probability maximizing] objective terms
time_limit = 10  # sec
noH = False

#Ophiopogon_japonicus scaffold 1
species = "Ophiopogon_japonicus"
scaffold = {
    "s0000009929": ["O=C1C(COc2ccccc12)Cc3ccc4OCOc4(c3)"]
}
scaffold_IDs = ["s0000000149"]
weights = [342.35, 370.36]
config_nums = [0] # the ith configuration recorded in the cIdx file (start from 0)

def load_scaffold_data(data_path):
    ''' 
    :param data_path: path of the target index file
    :return data object containing scaffold information
    '''

    lines = None
    with open(data_path) as f:
        lines = f.readlines()

    data = {}

    lptr = 0

    # read first line
    lptr += 1
    # read second line
    lptr += 1

    # read how many configurations
    cl = lines[lptr]
    data['num_config'] = int(cl.strip())
    lptr += 1

    # read and store the configurations
    configs = []
    for nconf in range(data['num_config']):
        cl = lines[lptr]
        cl_list = cl.split()

        num_spos = int(cl_list[0])
        sposs = [int(cl_list[ns]) for ns in range(1, num_spos+1)]

        configs.append((num_spos, sposs))
        lptr += 1

    data['configs'] = configs

    # read how many substituted positions (spos)
    cl = lines[lptr]
    data['num_spos'] = int(cl.strip())
    lptr += 1

    spos_sidechain_info = []

    for i in range(data['num_spos']):
        sp_sc_info = {}

        cl = lines[lptr]
        lptr += 1
        cl_list = cl.split()
        spos_id, spos_num_sidechain = int(cl_list[0]), int(cl_list[1])

        sp_sc_info['spos_id'] = spos_id

        sc_list = []
        for j in range(spos_num_sidechain):
            sc_info = {}
            cl = lines[lptr]
            lptr += 1
            cl_list = cl.split()
            sc_info['SMILES'] = cl_list[0]
            sc_info['freq'] = float(cl_list[1])
            sc_info['num_valence_electron'] = int(cl_list[2])
            sc_info['prob'] = float(cl_list[3])
            sc_info['MW'] = float(cl_list[4])
            sc_list.append(sc_info)

        sp_sc_info['sidechain_list'] = sc_list
        spos_sidechain_info.append(sp_sc_info)
    
    # sort based on spos id
    spos_sidechain_info.sort(key=lambda x: x['spos_id'])

    data['spos_sidechain_info'] = spos_sidechain_info

    return data

def prep_data(scaffold_ID):
    mol = Chem.MolFromSmiles(scaffold[scaffold_ID][0])
    molW = Descriptors.ExactMolWt(mol)
    
    data = load_scaffold_data(scaffold_ID+"_noH.cIdx")
    data['scaffold_weight'] = molW
    data['possible_molecular_weights'] = scaffold[scaffold_ID][1]
    return data

def m_create_bp(n, bp_matrix, bp_bias):
    term_bias = []
    if bp_bias != 0:
        term_bias = [{"c": bp_bias, "p": []}]
    terms_single = [{"c": bp_matrix[i, i], "p": [i]}
                    for i in range(n) if bp_matrix[i, i] != 0.0]
    terms_bonded = [{"c": bp_matrix[i, j], "p": [i, j]}
                    for i in range(n) for j in range(i+1, n) if bp_matrix[i, j] != 0.0]
    bp = {"terms": term_bias + terms_single + terms_bonded}
    return bp

def encode_dau_acceptable(tot_sc, num_sc, Wt, W, P, A=1, B=1, C=1):
    '''
    :param tot_sc: total number of possible sidechains
    :param num_sc: list of numbers of possible sidechains of each site
    :param Wt: total weight
    :param W: list of weights of possible sidechains of each site
    :param P: list of prob. of occurance of possible sidechains of each site
    :param A: weight for [one sidechain for one site] penalty terms
    :param B: weight for [weight equaling] penalty terms
    :param C: weight for [probability maximizing] objective terms
    :return: [binary_polynomial_object, penalty_binary_polynomial_object]
    '''

    W_flat = [wi for wc in W for wi in wc]
    P_flat = [pi for pc in P for pi in pc]
    log_P_flat = [math.log(pi, 2) for pi in P_flat]

    penalty_bp_matrix = np.zeros((tot_sc, tot_sc))
    penalty_bp_bias = 0.0
    # Actually bp_matrix can be 1D but we will fix bp_matrix as 2D
    # for now, so that it will be easier to combine penalty_bp_matrix
    # into bp_matrix if we want to!
    bp_matrix = np.zeros((tot_sc, tot_sc))
    bp_bias = 0.0

    # one sidechain for one site -> penalty
    cnt = 0
    for k in num_sc:
        for i in range(cnt, cnt+k):
            penalty_bp_matrix[i, i] -= 1 * A
            for j in range(i+1, cnt+k):
                penalty_bp_matrix[i, j] += 2 * A
        cnt += k
        penalty_bp_bias += 1 * A

    # weight equaling -> penalty
    for i in range(tot_sc):
        penalty_bp_matrix[i, i] += ((W_flat[i] ** 2) - 2 * Wt * W_flat[i]) * B
        for j in range(i+1, tot_sc):
            penalty_bp_matrix[i, j] += (2 * W_flat[i] * W_flat[j]) * B
    penalty_bp_bias += (Wt ** 2) * B

    # probability maximizing -> objective
    for i in range(tot_sc):
        # print(-log_P_flat[i] * C)
        bp_matrix[i, i] += -log_P_flat[i] * C

    binary_polynomial = m_create_bp(tot_sc, bp_matrix, bp_bias)
    penalty_binary_polynomial = m_create_bp(
        tot_sc, penalty_bp_matrix, penalty_bp_bias)

    return (binary_polynomial, penalty_binary_polynomial, num_sc)

def wrap(scaffold_ID, config_num, mw_num, A, B, C):
    data = prep_data(scaffold_ID)
    W = []  # list of weights of possible sidechains of each site
    P = []  # list of prob. of occurance of possible sidechains of each site
    num_sc = []  # list of numbers of possible sidechains of each site
    for spos in data['spos_sidechain_info']:
        spos_id = spos['spos_id']
        sidechain_list = spos['sidechain_list']
    
        tmp_w = []
        tmp_p = []
        for sidechain in sidechain_list:
    
            prob = sidechain['prob']
    
            tmp_w.append(sidechain['MW'])
            tmp_p.append(sidechain['prob'])
    
        W.append(tmp_w)
        P.append(tmp_p)
        num_sc.append(len(tmp_w))
    
    config = data['configs'][config_num]
    
    num_spos = config[0]
    sposs = config[1]
    
    c_W = list(W[i] for i in sposs)  # [[w11, w12, w13], [w21, w22], ...]
    c_P = list(P[i] for i in sposs)
    c_num_sc = list(num_sc[i] for i in sposs)  # [3, 2, ...]
    c_tot_sc = sum(c_num_sc)
    
    H_weight = 1.007825032
    Wt = data['possible_molecular_weights'][mw_num] - (data['scaffold_weight'] - num_spos * H_weight)
    
    return encode_dau_acceptable(c_tot_sc, c_num_sc, Wt, c_W, c_P, A, B, C)

for scaffold_ID in scaffold_IDs:
    for weight in weights:
        for config_num in config_nums:  
            ((bp, penalty_bp, num_sc, W_sidechains), data) = wrap(scaffold_ID, weight, config_num, A, B, C, noH)
            
            for i in range(1):
                data_configuration = data_config(data, config_num)
                
                job_id = wrap_DAU(bp, penalty_bp, num_sc, API_Key, API_Access_URL, headers, proxies, time=time_limit)
                print(job_id)
                time.sleep(time_limit+40)
                
                response = requests.get(API_Access_URL + '/da/v3/async/jobs/result/' + job_id, headers=headers, proxies=proxies)
                print(response.json())
                out = interprete_response(response.json(), data_configuration, W_sidechains)
                delete = requests.delete(API_Access_URL + '/da/v3/async/jobs/result/' + job_id, headers=headers, proxies=proxies)
                
                if len(out) > 0:
                    filename = f"{species}_{scaffold_ID}_config_{config_num}_peak_{weight}_withH_{i}"
                    
                    write_csv("results/" + filename, data_configuration, out)
                    print("weight:", weight, "config:", config_num, "Yes")
                else:
                    print("weight:", weight, "config:", config_num, "No")
