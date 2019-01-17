function valid_idx= gate_meas_smc(z,P_G,model,x,w,ngrid)

valid_idx = [];
zlength = size(z,2); if zlength==0, z_gate= []; return; end
plength = size(x,2);

%edges for cells
edges= cell(model.z_dim,1);
for d=1:model.z_dim
    edges{d}= linspace(model.range_c(d,1),model.range_c(d,2),ngrid+1);
end

%the histogram for predicted track
eta= gen_observation_fn(model,x,'noiseless');
[wphist xphist hpedges hpmid xploc]= histcn(w,eta',edges);
[xxxnmlp,xxxidxp]= sort(wphist(:)/sum(wphist(:)),'descend');
[yyyidxp]= find(cumsum(xxxnmlp)>P_G,1);
gateidx=xxxidxp(1:yyyidxp);

%the gating for predicted track
[wzhist zfhist hzedges hzmid zfloc]= histcn(ones(zlength,1),z',edges);
for emm=1:zlength
    zedemmloc=cell(model.z_dim,1);
    for d=1:model.z_dim
        zedemmloc{d}= zfloc(emm,d);
    end
    if all(zfloc(emm,:)>0)
        if any(sub2ind([ngrid ngrid],zedemmloc{:})==gateidx)
            valid_idx= unique_faster([valid_idx emm]);
        end
    end
end
valid_idx=valid_idx(:)';