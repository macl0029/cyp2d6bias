import delimited "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data_pm.csv", varnames(1) clear 
drop if effect=="NA"

g excl=0
replace excl=1 if name=="Thompson" | name=="Wegman" | name=="Chamnanphon"|name=="Teh"

replace name="Goetz" if name=="Goetzb"
replace author="Goetz, 2013" if author=="Goetzb, 2013"


destring effect,replace force
destring up,replace force
destring low,replace force

destring effc,replace force
destring upc,replace force
destring lowc,replace force

destring effloh,replace force
destring uploh,replace force
destring lowloh,replace force

destring yearpu, replace force
g yeard=substr(datayearrange,1,4)
destring yeard, replace force
destring age, replace force
destring post, replace force
destring outcome_cat, replace force

replace followup="10" if followup=="<10"
replace followup="12" if followup=="<12y"
replace followup="14" if followup=="<14"
replace followup="3" if followup=="<3"
destring followup, replace force

g cohortstudy=studytype=="cohort"
g imbias=immortaltime=="1"
recode imbias .=0


g source=(dnasource=="Tumor")

g leff=log(effect)
g lup=log(up)
g llow=log(low)
g leff_c=log(effc)
g lup_c=log(upc)
g llow_c=log(lowc)
g leff_l=log(effloh)
g llow_l=log(lowloh)
g lup_l=log(uploh)

*Crude results
meta set leff_c llow_c lup_c, civartolerance(.5) studylabel(author) 
meta summarize, eform cformat(%9.2f)

*LOH only
meta set leff_l llow_l lup_l, civartolerance(.5) studylabel(author)
meta summarize , eform cformat(%9.2f)

*loh and inc
meta set leff llow lup ,civartolerance(.5) studylabel(author) 
meta summarize , eform cformat(%9.2f)
meta forestplot , esrefline nullrefline sort(yearp) eform(Risk Ratio) noosigtest noohet noohom 

*eliminate studies with impossible values
meta summarize if excl==0, eform cformat(%9.2f)
meta forestplot if excl==0, esrefline nullrefline sort(yearp) eform(Risk Ratio) noosigtest noohet noohom 

*loh and inc - only if have genotype/phenotype data
g nodata=1 if name=="Markkula"|name=="Thompson"|name=="Wegman"|name=="Goetz"|name=="Gor"|name=="Rae"

meta set leff llow lup ,civartolerance(.9) studylabel(author) 
meta summarize if nodata!=1, eform cformat(%9.2f)
	meta forestplot if nodata!=1 , esrefline nullrefline sort(yearp) eform(Risk Ratio) noosigtest noohet noohom


*Publication bias tests
meta bias, egger
meta bias, begg

meta regress _cons , eform
meta regress source, eform
meta regress age , eform
meta regress yearpu , eform
meta regress yeard , eform
meta regress post , eform
meta regress followup , eform
meta regress cohortstudy , eform
meta regress i.outcome_cat ,eform
meta summarize , subgroup(outcome_cat) eform

meta summarize, eform subgroup(source)

****Repeat for IM vs NM - eur
import delimited "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data_im.csv", varnames(1) clear 
drop if effect=="NA"
replace name="Goetz" if name=="Goetzb"
replace author="Goetz, 2013" if author=="Goetzb, 2013"


g excl=0
replace excl=1 if name=="Thompson" | name=="Wegman" | name=="Chamnanphon"|name=="Teh"

destring effect,replace force
destring up,replace force
destring low,replace force

destring effc,replace force
destring upc,replace force
destring lowc,replace force

destring effloh,replace force
destring uploh,replace force
destring lowloh,replace force

destring yearpu, replace force
g yeard=substr(datayearrange,1,4)
destring yeard, replace force
destring age, replace force
destring post, replace force
destring outcome_cat, replace force

replace followup="10" if followup=="<10"
replace followup="12" if followup=="<12y"
replace followup="14" if followup=="<14"
replace followup="3" if followup=="<3"
destring followup, replace force

g cohortstudy=studytype=="cohort"
g imbias=immortaltime=="1"
recode imbias .=0

g source=(dnasource=="Tumor")


g leff=log(effect)
g lup=log(up)
g llow=log(low)
g leff_c=log(effc)
g lup_c=log(upc)
g llow_c=log(lowc)
g leff_l=log(effloh)
g llow_l=log(lowloh)
g lup_l=log(uploh)

*Crude results
meta set leff_c llow_c lup_c, civartolerance(.3) studylabel(author)
meta summarize , eform cformat(%9.2f)

*LOH only
meta set leff_l llow_l lup_l, civartolerance(.3) studylabel(author)
meta summarize , eform cformat(%9.2f)

*loh and inc
meta set leff llow lup ,civartolerance(.9) studylabel(author) 
meta summarize , eform cformat(%9.2f)
meta forestplot , esrefline nullrefline sort(yearp) eform(Risk Ratio) noosigtest noohet noohom



	



meta bias, egger
meta bias, begg

meta regress _cons , eform
meta regress source, eform
meta regress age , eform
meta regress yearpu , eform
meta regress yeard , eform
meta regress post , eform
meta regress followup , eform
meta regress cohortstudy , eform
meta regress i.outcome_cat ,eform
meta summarize , subgroup(outcome_cat) eform



****Repeat for IM vs NM - as
import delimited "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data_imas.csv", varnames(1) clear 
drop if effect=="NA"




destring effect,replace force
destring up,replace force
destring low,replace force

destring effc,replace force
destring upc,replace force
destring lowc,replace force

destring effloh,replace force
destring uploh,replace force
destring lowloh,replace force

destring yearpu, replace force
g yeard=substr(datayearrange,1,4)
destring yeard, replace force
destring age, replace force
destring post, replace force
destring outcome_cat, replace force
destring followup, replace force

g cohortstudy=studytype=="cohort"
g imbias=immortaltime=="1"
recode imbias .=0


g source=(dnasource=="Tumor")


g leff=log(effect)
g lup=log(up)
g llow=log(low)
g leff_c=log(effc)
g lup_c=log(upc)
g llow_c=log(lowc)
g leff_l=log(effloh)
g llow_l=log(lowloh)
g lup_l=log(uploh)

*Crude results
meta set leff_c llow_c lup_c, civartolerance(.3) studylabel(author)
meta summarize, eform cformat(%9.2f)

*LOH only
meta set leff_l llow_l lup_l, civartolerance(.3) studylabel(author)
meta summarize , eform cformat(%9.2f)

*loh and inc
meta set leff llow lup ,civartolerance(.9) studylabel(author) 
meta summarize , eform cformat(%9.2f)
meta forestplot , esrefline nullrefline sort(yearp) eform(Risk Ratio) noosigtest noohet noohom


meta bias, egger
meta bias, begg

meta regress _cons , eform
meta regress source, eform
meta regress age , eform
meta regress yearpu , eform
meta regress yeard , eform
meta regress post , eform
meta regress followup , eform
meta regress cohortstudy , eform
meta regress i.outcome_cat ,eform
meta summarize , subgroup(outcome_cat) eform



****Repeat for PMIM vs NM 

import delimited "~/Dropbox/Projects - Dropbox/reproducibility/mine/cyp2d6 misclass/meta_data_pmim.csv", varnames(1) clear 
drop if effect=="NA"

replace name="Goetz" if name=="Goetzb"
replace author="Goetz, 2013" if author=="Goetzb, 2013"



destring effect,replace force
destring up,replace force
destring low,replace force

destring effc,replace force
destring upc,replace force
destring lowc,replace force

destring effloh,replace force
destring uploh,replace force
destring lowloh,replace force

destring yearpu, replace force
g yeard=substr(datayearrange,1,4)
destring yeard, replace force
destring age, replace force
destring post, replace force
destring outcome_cat, replace force
destring followup, replace force

g cohortstudy=studytype=="cohort"
g imbias=immortaltime=="1"
recode imbias .=0



g leff=log(effect)
g lup=log(up)
g llow=log(low)
g leff_c=log(effc)
g lup_c=log(upc)
g llow_c=log(lowc)
g leff_l=log(effloh)
g llow_l=log(lowloh)
g lup_l=log(uploh)

*Crude results
meta set leff_c llow_c lup_c, civartolerance(.3) studylabel(author)
meta summarize, eform cformat(%9.2f)

*LOH only
meta set leff_l llow_l lup_l, civartolerance(.3) studylabel(author)
meta summarize , eform cformat(%9.2f)

*loh and inc
meta set leff llow lup ,civartolerance(.9) studylabel(author) 
meta summarize , eform cformat(%9.2f)
meta forestplot , esrefline nullrefline sort(yearp) eform(Risk Ratio) noosigtest noohet noohom


meta bias, egger
meta bias, begg

meta regress _cons , eform
meta regress age , eform
meta regress yearpu , eform
meta regress yeard , eform
meta regress post , eform
meta regress followup , eform
meta regress cohortstudy , eform
meta regress i.outcometype ,eform
meta summarize , subgroup(outcome_cat) eform



