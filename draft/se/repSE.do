* Script to make Extended Data Table 1

clear all
// set mem 1G
// set matsize 10000
// set maxvar 10000
cd "/Users/tombearpark/Library/CloudStorage/Dropbox/BP_2023_fesearch/data/BurkeHsiangMiguel2015_Replication"

// postfile robust mod b1 b2 using data/output/Coefficients_robustness, replace
use data/input/GrowthClimateDataset, clear
xtset iso_id year
gen temp = UDel_temp_popweight
gen temp2 = temp*temp
gen precip = UDel_precip_popweight/1000 //rescaling so precip coeffs woudl be legible
gen precip2 = precip*precip
cap est drop reg*


reg growthWDI temp temp2 precip precip2 i.year _yi_* _y2_* i.iso_id , cluster(iso_id)  //full sample
est sto m1

reg growthWDI temp temp2 precip precip2 i.year _yi_* _y2_* i.iso_id , cluster(iso_id, hc3)  //full sample
est sto m2


esttab m1 m2, keep(temp temp2) se


