% 
hinfo_bl = hdf5info('20190927_A42_Ch1_S1_Ch2_M1_baseline_0001');
dset_bl = hdf5read(hinfo_bl.GroupHierarchy.Groups(2).Datasets(1));
dset_bl_M1 = dset_bl(:,2);
filta1 = filter(alpha,1,dset2);
filtb1 = filter(beta,1,dset2);
filtc1 = filter(gamma,1,dset2);
filtt1 = filter(theta,1,dset2);

hinfo_pilo = hdf5info('20190927_A42_Ch1_S1_Ch2_M1_pilo-2_0001');
dset_pilo = hdf5read(hinfo_pilo.GroupHierarchy.Groups(2).Datasets(1));
dset_pilo_M1 = dset_pilo(:,2);
