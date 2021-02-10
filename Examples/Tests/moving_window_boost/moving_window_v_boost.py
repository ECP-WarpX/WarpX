import os
import yt

dlist0 = os.listdir('./diags')
dlist = sorted(dlist0)

kk = 0
fn = './diags/'+dlist[kk]
print(fn)
ds=yt.load(fn)
t1=ds.current_time.v[()]
z1max=ds.domain_right_edge.v[2]

kk = 1 
fn = './diags/'+dlist[kk]
print(fn)
ds=yt.load(fn)
t2=ds.current_time.v[()]
z2max=ds.domain_right_edge.v[2]

vv=(z2max-z1max)/(t2-t1)/299792458.
print('moving_window_v in boosted frame =', vv)
