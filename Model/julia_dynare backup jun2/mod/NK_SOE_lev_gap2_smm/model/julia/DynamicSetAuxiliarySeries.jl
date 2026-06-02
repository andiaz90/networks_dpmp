function dynamic_set_auxiliary_series!(ds, params)
#
# Computes auxiliary variables of the dynamic model
#
@inbounds begin
ds.AUX_ENDO_LAG_20_1 .=lag(ds.om_g);
ds.AUX_ENDO_LAG_21_1 .=lag(ds.om_s);
ds.AUX_ENDO_LAG_138_1 .=lag(ds.A_1);
ds.AUX_ENDO_LAG_139_1 .=lag(ds.A_2);
ds.AUX_ENDO_LAG_140_1 .=lag(ds.A_3);
ds.AUX_ENDO_LAG_141_1 .=lag(ds.A_4);
ds.AUX_ENDO_LAG_142_1 .=lag(ds.A_5);
ds.AUX_ENDO_LAG_143_1 .=lag(ds.A_6);
ds.AUX_ENDO_LAG_144_1 .=lag(ds.A_7);
ds.AUX_ENDO_LAG_145_1 .=lag(ds.A_8);
ds.AUX_ENDO_LAG_146_1 .=lag(ds.A_9);
ds.AUX_ENDO_LAG_147_1 .=lag(ds.A_10);
ds.AUX_ENDO_LAG_148_1 .=lag(ds.A_11);
ds.AUX_ENDO_LAG_149_1 .=lag(ds.A_12);
end
end
