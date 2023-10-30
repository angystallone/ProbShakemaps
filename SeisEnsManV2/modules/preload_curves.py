import os
import sys
import glob
import h5py
from pymatreader import read_mat
import numpy as np
from preload import check_if_path_exists

def reallocate_curves(**kwargs):

    Config    = kwargs.get('cfg', None)
    args      = kwargs.get('args', None)
    c_files   = kwargs.get('curve_files', None)
    name      = kwargs.get('name', None)

    #curves_py_folder  =  Config.get('pyptf',  'curves')
    #curves_py_folder  =  Config.get('pyptf',  'curves_gl_16')

    curves_py_folder  =  Config.get('pyptf',  'h_curves')

    list_out = []
    for i in range(0,int(Config.get('ScenariosList', 'nr_regions'))):
        d = "%03d" % (i+1)
        def_name = curves_py_folder + os.sep + name + d + '-empty.hdf5'
        #def_name = curves_py_folder + os.sep + name + d + '-empty.npy'
        list_out.append(def_name)

    for i in range(len(c_files)):
        ref_nr = int(c_files[i].split(name)[-1][0:3])
        list_out[ref_nr-1] = c_files[i]

    return list_out

def load_hazard_values(**kwargs):

    Config    = kwargs.get('cfg', None)
    args      = kwargs.get('args', None)
    in_memory = kwargs.get('in_memory', False)

    """
    Preload from mat file very time consuming, about 10 hours, so args.preload disabled
    for this section. Only args.preload_scenarios
    """



    # inizialize empty dict, containing the list of the ps and bs scenarios
    scenarios    = dict()

    # Load scenarios path
    # for hdf5
    #curves_py_folder  =  Config.get('pyptf',  'curves')
    # for npy
    #curves_py_folder  =  Config.get('pyptf',  'curves_gl_16')
    curves_py_folder  =  Config.get('pyptf',  'h_curves')
    #curves_mat_folder =  Config.get('tsumaps','curves')


    # Load ps and bs in pypath
    py_os_ps_curves = sorted(glob.glob(curves_py_folder + os.sep + Config.get('pyptf','os_ps_curves_file_names') + '*'))
    py_os_bs_curves = sorted(glob.glob(curves_py_folder + os.sep + Config.get('pyptf','os_bs_curves_file_names') + '*'))
    py_af_ps_curves = sorted(glob.glob(curves_py_folder + os.sep + Config.get('pyptf','af_ps_curves_file_names') + '*'))
    py_af_bs_curves = sorted(glob.glob(curves_py_folder + os.sep + Config.get('pyptf','af_bs_curves_file_names') + '*'))
    py_gl_ps_curves = sorted(glob.glob(curves_py_folder + os.sep + Config.get('pyptf','gl_ps_curves_file_names') + '*'))
    py_gl_bs_curves = sorted(glob.glob(curves_py_folder + os.sep + Config.get('pyptf','gl_bs_curves_file_names') + '*'))

    py_os_ps_curves = reallocate_curves(curve_files=py_os_ps_curves, args=args, cfg=Config, name=Config.get('pyptf','os_ps_curves_file_names'))
    py_os_bs_curves = reallocate_curves(curve_files=py_os_bs_curves, args=args, cfg=Config, name=Config.get('pyptf','os_bs_curves_file_names'))
    py_af_ps_curves = reallocate_curves(curve_files=py_af_ps_curves, args=args, cfg=Config, name=Config.get('pyptf','af_ps_curves_file_names'))
    py_af_bs_curves = reallocate_curves(curve_files=py_af_bs_curves, args=args, cfg=Config, name=Config.get('pyptf','af_bs_curves_file_names'))
    py_gl_ps_curves = reallocate_curves(curve_files=py_gl_ps_curves, args=args, cfg=Config, name=Config.get('pyptf','gl_ps_curves_file_names'))
    py_gl_bs_curves = reallocate_curves(curve_files=py_gl_bs_curves, args=args, cfg=Config, name=Config.get('pyptf','gl_bs_curves_file_names'))


#    if (len(py_os_ps_curves) == 0 or args.preload_curves == 'Yes'):
#
#        path_py_exist = check_if_path_exists(path=curves_py_folder, create=True)
#        print("PreLoad ps-os Curves: npy conversion     <------ ", curves_mat_folder)
#        mat_os_ps_curves = sorted(glob.glob(curves_mat_folder + os.sep + Config.get('tsumaps','os_ps_curves_file_names') + '*'))
#        py_os_ps_curves  = mat_curves_to_py_curves(mat_path = curves_mat_folder,
#                                                   py_path  = curves_py_folder,
#                                                   files    = mat_os_ps_curves,
#                                                  mat_key  = 'osVal_PS',
#                                                   cfg      = Config)
#
#    if (len(py_os_bs_curves) == 0 or args.preload_curves == 'Yes'):
#
#        path_py_exist = check_if_path_exists(path=curves_py_folder, create=True)
#        print("PreLoad ps-os Curves: npy conversion     <------ ", curves_mat_folder)
#        mat_os_bs_curves = sorted(glob.glob(curves_mat_folder + os.sep + Config.get('tsumaps','os_bs_curves_file_names') + '*'))
#        py_os_bs_curves  = mat_curves_to_py_curves(mat_path = curves_mat_folder,
#                                                   py_path  = curves_py_folder,
#                                                   files    = mat_os_bs_curves,
#                                                   mat_key  = 'osVal_BS',
#                                                   cfg      = Config)
#
#    if (len(py_af_ps_curves) == 0 or args.preload_curves == 'Yes'):
#
#        path_py_exist = check_if_path_exists(path=curves_py_folder, create=True)
#        print("PreLoad ps-af Curves: npy conversion     <------ ", curves_mat_folder)
#        mat_af_ps_curves = sorted(glob.glob(curves_mat_folder + os.sep + Config.get('tsumaps','af_ps_curves_file_names') + '*'))
#        py_af_ps_curves  = mat_curves_to_py_curves(mat_path = curves_mat_folder,
#                                                   py_path  = curves_py_folder,
#                                                   files    = mat_af_ps_curves,
#                                                   mat_key  = 'afVal_PS',
#                                                   cfg      = Config)
#
#    if (len(py_af_bs_curves) == 0 or args.preload_curves == 'Yes'):
#
#        path_py_exist = check_if_path_exists(path=curves_py_folder, create=True)
#        print("PreLoad bs-af Curves: npy conversion     <------ ", curves_mat_folder)
#        mat_af_bs_curves = sorted(glob.glob(curves_mat_folder + os.sep + Config.get('tsumaps','af_bs_curves_file_names') + '*'))
#        py_af_bs_curves  = mat_curves_to_py_curves(mat_path = curves_mat_folder,
#                                                   py_path  = curves_py_folder,
#                                                   files    = mat_af_bs_curves,
#                                                   mat_key  = 'afVal_BS',
#                                                   cfg      = Config)
#
#
#    ##########
#    if (len(py_gl_ps_curves) == 0 or args.preload_curves == 'Yes'):
#
#        path_py_exist = check_if_path_exists(path=curves_py_folder, create=True)
#        print("PreLoad ps-gl Curves: npy conversion     <------ ", curves_mat_folder)
#        mat_gl_ps_curves = sorted(glob.glob(curves_mat_folder + os.sep + Config.get('tsumaps','gl_ps_curves_file_names') + '*'))
#        py_gl_ps_curves  = mat_curves_to_py_curves(mat_path = curves_mat_folder,
#                                                   py_path  = curves_py_folder,
#                                                   files    = mat_gl_ps_curves,
#                                                   mat_key  = 'glVal_PS',
#                                                   cfg      = Config)
#
#    if (len(py_gl_bs_curves) == 0 or args.preload_curves == 'Yes'):
#
#        path_py_exist = check_if_path_exists(path=curves_py_folder, create=True)
#        print("PreLoad ps-gl Curves: npy conversion     <------ ", curves_mat_folder)
#        mat_gl_bs_curves = sorted(glob.glob(curves_mat_folder + os.sep + Config.get('tsumaps','gl_bs_curves_file_names') + '*'))
#        py_gl_bs_curves  = mat_curves_to_py_curves(mat_path = curves_mat_folder,
#                                                   py_path  = curves_py_folder,
#                                                   files    = mat_gl_bs_curves,
#                                                   mat_key  = 'glVal_BS',
#                                                   cfg      = Config)

    hazard_curves_files = dict()
    hazard_curves_files['gl_ps'] = py_gl_ps_curves
    hazard_curves_files['gl_bs'] = py_gl_bs_curves
    hazard_curves_files['os_ps'] = py_os_ps_curves
    hazard_curves_files['os_bs'] = py_os_bs_curves
    hazard_curves_files['af_ps'] = py_af_ps_curves
    hazard_curves_files['af_bs'] = py_af_bs_curves


    return hazard_curves_files

def mat_curves_to_py_curves(**kwargs):

    curves_py_folder  = kwargs.get('py_path', None)
    curves_mat_folder = kwargs.get('mat_path', None)
    files             = kwargs.get('files', None)
    vmat              = kwargs.get('vmat', None)
    mat_key           = kwargs.get('mat_key', None)

    npy_curves_files = []

    for i in range(len(files)):

        npy_file = files[i].replace(curves_mat_folder,curves_py_folder).replace('.mat','.hdf5')

        #if (os.path.isfile(npy_file)):
        #    continue
        print('   ', npy_file, ' <--- ', files[i])
        py_dict  = read_mat(files[i])
        hf = h5py.File(npy_file, 'w')
        hf.create_dataset(mat_key, data=py_dict[mat_key])
        hf.close()
        #potrebbe aiutare?
        #del(hf)

        #"""
        #key      = [*py_dict][0]
        #np.save(npy_file, py_dict[mat_key]) #, allow_pickle=True)
        #"""
        npy_curves_files.append(npy_file)

        #print(py_dict['osVal_PS'][0:3][:])
        #print(type(py_dict['osVal_PS']))
        #sys.exit()

        """
        try:
            py_dict = read_mat(files[i])
        except:
            continue

        try:
            #np.save(npy_file, py_dict, allow_pickle=True)
            done = save_dict(npy=npy_file, dict=py_dict, cfg=Config)
        except:
            continue

        npy_scenarios_files.append(npy_file)
        print('  --> OK: ', npy_file)
        """
    #print(npy_curves_files)


    return npy_curves_files
