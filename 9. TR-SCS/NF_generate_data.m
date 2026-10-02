function NF = NF_generate(LORs_frac, crystal_lut, PET_geom, NF_ax_data, NF_tr_data)
% Modified to accept data structures instead of filenames

% Extract data from structures
b_ax = NF_ax_data.b_ax;
g_ax = NF_ax_data.g_ax;
% Check where epsilon is located (usually in one of them)
if isfield(NF_ax_data, 'epsilon')
    epsilon = NF_ax_data.epsilon;
else
    epsilon = NF_tr_data.epsilon;
end

g_tr = NF_tr_data.g_tr;

L = PET_geom(:,1).*PET_geom(:,3).*PET_geom(:,4);
M = PET_geom(:,2).*PET_geom(:,5);

% Pass crystal_lut directly to avoid repeated file reading
lor = event2lor_data(LORs_frac, crystal_lut, PET_geom);

% module Crystal ID logic
histo_add(:,1) = mod(lor(:,3),PET_geom(1,1)*PET_geom(1,3)).*(lor(:,1)<(L(1)*M(1)));
histo_add(:,2) = mod(lor(:,6),PET_geom(1,1)*PET_geom(1,3)).*(lor(:,4)<(L(1)*M(1)));

lor(:,7)= (((lor(:,6)-lor(:,3)) < 0 )*sum(L)  + lor(:,6)-lor(:,3))+1;

NF_all=zeros(length(lor),6);
% [b_ax(u)  b_ax(v)  g_ax(uv)  g_tr(r)  f_tr(r)  e(ui)  e(vj)]
NF_all(:,1) = b_ax(lor(:,2)+1);
NF_all(:,2) = b_ax(lor(:,5)+1);
NF_all(:,3) = g_ax(sub2ind(size(g_ax),lor(:,2)+1,lor(:,5)+1));
NF_all(:,4) = g_tr(sub2ind(size(g_tr),lor(:,7),histo_add(:,1)+1));
NF_all(:,5) = epsilon(sub2ind(size(epsilon),lor(:,2)+1,lor(:,3)+1));
NF_all(:,6) = epsilon(sub2ind(size(epsilon),lor(:,5)+1,lor(:,6)+1));
NF=prod(NF_all,2);

end
