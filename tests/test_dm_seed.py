import numpy as np
import pytest
from phasebuild.dm_seed import canonicalize

def test_one_ulp_variants_share_map_seed_without_mutating_originals():
 phase=np.array([-180.,-72.24123,0.,137.25235,180.],np.float32);fom=np.array([0.,.23156722,.5,.8806152,1.],np.float32)
 original=phase.copy();changed=phase.copy();changed[1]=np.nextafter(changed[1],np.float32(0));other=fom.copy();other[1]=np.nextafter(other[1],np.float32(1))
 p,w,a=canonicalize(phase,fom);q,v,b=canonicalize(changed,other)
 np.testing.assert_array_equal(p,q);np.testing.assert_array_equal(w,v);np.testing.assert_array_equal(phase,original)
 assert a['max_phase_change_degrees']<=.00055 and a['max_fom_change']<=.00000055
 assert p[0]==-180 and p[-1]==180 and w[0]==0 and w[-1]==1

@pytest.mark.parametrize('phase,fom',[([0],[1.01]),([float('nan')],[.5]),([0],[float('inf')]),([[0]],[[.5]]),([0,1],[.5])])
def test_bad_seed_rejected(phase,fom):
 with pytest.raises(ValueError,match='DM seed'):canonicalize(phase,fom)
