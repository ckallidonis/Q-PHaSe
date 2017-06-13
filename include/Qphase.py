#
# include/Qphase.py
#
# Various global definitions
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#----------------------------------------

# Flavor combinations and correlation function direction list
flav_list2pt = ['ppm','pmm']
drct_list2pt = ['fwd','bwd']
flav_list3pt = ['up','dn']
flcombs_conn = ['up','dn','IS','IV']

# Supported Gamma-matrix bases
Gamma_bases = {'QuaHoG':'gammas_QuaHoG' , 'QUDA_UKQCD':'gammas_QUDA_UKQCD', 'QUDA_DeGrand':'gammas_QUDA_DeGrand_Rossi'}

# Supported operators
supported_operators = ('axial-vec' , 'scalar')

# Projector tags for each basis ( Convention: (x,y,xyz) )
projector_tags3pt = { 'axial-vec': {'QuaHoG_polarized':('P4','P5','P3')} , 'scalar': {'QuaHoG_unpolarized':('P0',)} }

# Insertion current tags for each operator and basis ( Convention: (x,y,z) ) *** Assume that the three-point function is in the twisted basis in QuaHoG basis!!! ***
operator_tags3pt = { 'axial-vec': {'QuaHoG':('=loc:g5gx=','=loc:g5gy=','=loc:g5gz=')} , 'scalar': {'QuaHoG':('=loc:g5=',)} }

# This helps flip the sign for the scalar and the tensor charges if we are in the twisted basis
optr_sign3pt = {'QuaHoG': {'axial-vec': +1.0 , 'scalar': -1.0}  }

#---------------------------------------------------------

def check_support(basis_tag,optr_tag):
    
    if( basis_tag != 'QuaHoG'): # Only QuaHoG basis supported for now
        print('get_thrp: Only QuaHoG basis supported for now! Exiting.')
        sys.exit()
        
    if( optr_tag not in supported_operators):
        print('Operator %s not supported! Exiting.' % (optr_tag))
        sys.exit()
#---------------------------------------------------------

    
def get_3ptproj_tag(basis_tag,optr_tag):

    proj_pol = ''
    if(optr_tag == 'axial-vec'):
        proj_pol = 'polarized'
    elif(optr_tag == 'scalar'):
        proj_pol = 'unpolarized'
    else:
        print( 'get_3ptproj_tag: Please add the operator %s inside get_3ptproj_tag!' % (optr_tag) ) 
        
    return basis_tag + '_' + proj_pol
#---------------------------------------------------------
