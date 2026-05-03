import glob,csv,math,json,os

G=6.67430e-11

def read_snap(path):
    rows=[]
    with open(path) as f:
        r=csv.DictReader(f)
        for row in r:
            rows.append({k:float(v) for k,v in row.items() if k in ['x','y','vx','vy','mass']})
    return rows

def accel(state,bh_mass,soft,star_star):
    n=len(state)
    ax=[0.0]*n; ay=[0.0]*n
    eps2=soft*soft
    for i,a in enumerate(state):
        dx=-a['x']; dy=-a['y']; r2=dx*dx+dy*dy+eps2; inv=r2**-1.5
        ax[i]+=G*bh_mass*dx*inv; ay[i]+=G*bh_mass*dy*inv
    if star_star:
        for i,a in enumerate(state):
            for j,b in enumerate(state):
                if i==j: continue
                dx=b['x']-a['x']; dy=b['y']-a['y']; r2=dx*dx+dy*dy+eps2; inv=r2**-1.5
                ax[i]+=G*b['mass']*dx*inv; ay[i]+=G*b['mass']*dy*inv
    return ax,ay

def case(name,n_dir,t_dir,bh,soft,star_star):
    nfiles=sorted(glob.glob(f"{n_dir}/snapshot_*.csv"))
    tfiles=sorted(glob.glob(f"{t_dir}/snapshot_*.csv"))
    m=min(len(nfiles),len(tfiles))
    max_ad=max_ar=max_pd=max_vd=0.0
    for k in range(m):
        ns=read_snap(nfiles[k]); ts=read_snap(tfiles[k])
        nax,nay=accel(ns,bh,soft,star_star); tax,tay=accel(ts,bh,soft,star_star)
        for i in range(min(len(ns),len(ts))):
            dax=tax[i]-nax[i]; day=tay[i]-nay[i]; d=math.hypot(dax,day)
            nmag=math.hypot(nax[i],nay[i])
            max_ad=max(max_ad,d)
            if nmag>0: max_ar=max(max_ar,d/nmag)
            dp=math.hypot(ts[i]['x']-ns[i]['x'],ts[i]['y']-ns[i]['y']); max_pd=max(max_pd,dp)
            dv=math.hypot(ts[i]['vx']-ns[i]['vx'],ts[i]['vy']-ns[i]['vy']); max_vd=max(max_vd,dv)
    return dict(case=name,snapshots_compared=m,max_abs_accel_diff=max_ad,max_rel_accel_diff=max_ar,max_pos_diff=max_pd,max_vel_diff=max_vd)

results=[
case('case1_single_source','outputs/val_case1_newt','outputs/val_case1_tpf',1.98847e36,0.0,False),
case('case2_galaxy_starstar_false','outputs/val_case2_newt','outputs/val_case2_tpf',1.98847e36,1e15,False),
case('case3_smallN_starstar_true','outputs/val_case3_newt','outputs/val_case3_tpf',1.98847e35,5e14,True),
]
os.makedirs('validate/geodesic_correspondence_newtonian/artifacts',exist_ok=True)
with open('validate/geodesic_correspondence_newtonian/artifacts/metrics.json','w') as f: json.dump(results,f,indent=2)
print(json.dumps(results,indent=2))
