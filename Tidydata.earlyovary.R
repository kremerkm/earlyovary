library(tidyverse)
library(SEER2R)
library(zoo)
library(survival)
library(survminer)
library(kableExtra)
library(RColorBrewer)

#Read SEER case files for ovarian cancer from 2004-2016 into data.frame
ovary.df <- read.SeerStat("ovcase_2004-2016.dic", UseVarLabelsInData = TRUE)

#Change data.frame to tibble
ovary.tib <- as_tibble(ovary.df)

#Rename columns
ovary.tib = rename(ovary.tib,
            Race = "Race_and_origin_recode_NHW_NHB_NHAIAN_NHAPI_Hispanic",
            Histo = "Histologic_Type_ICDO3",
            Lat = "Laterality",
            Summ_Stage = "SEER_Combined_Summary_Stage_2000_2004",
            Stage_Group = "Derived_AJCC_Stage_Group_6th_ed_20042015",
            T_Stage = "Derived_AJCC_T_6th_ed_20042015",
            N_Stage = "Derived_AJCC_N_6th_ed_20042015",
            M_Stage = "Derived_AJCC_M_6th_ed_20042015",
            Nodes_Ex = "Regional_nodes_examined_1988",
            Nodes_Pos = "Regional_nodes_positive_1988",
            Chemo = "Chemotherapy_recode_yes_no/unk",
            DiagYear = "Year_of_diagnosis",
            SurvMonths = "Survival_months",
            COD = "SEER_causespecific_death_classification",
            DiagAge = "Age_at_diagnosis",
            CutoffStatus = "Vital_status_recode_study_cutoff_used",
            DiagMonth = "Month_of_diagnosis_recode",
            CutoffStatus2 = "End_Calc_Vital_Status_Adjusted",
            MonthsFollowed = "Number_of_Intervals_Calculated")

#Remove columns
ovary.tib = select(ovary.tib, -c(Summary_stage_2000_1998, Sex, County, CutoffStatus2, Begin_Calc_Age_Adjusted))

#Change Grades to numeric values
ovary.tib$Grade <- ifelse(ovary.tib$Grade == "Poorly differentiated; Grade III", 3,
        ifelse(ovary.tib$Grade == "Moderately differentiated; Grade II", 2, 
        ifelse(ovary.tib$Grade == "Well differentiated; Grade I", 1,
        ifelse(ovary.tib$Grade == "Undifferentiated; anaplastic; Grade IV", 4, "Unk"))))

#Change Lat to Right/Left/Bil
ovary.tib$Lat <- ifelse(ovary.tib$Lat == "Right - origin of primary", "R",
        ifelse(ovary.tib$Lat == "Left - origin of primary", "L", 
        ifelse(ovary.tib$Lat == "Bilateral, single primary", "Bil", "Unk")))

#Change Race to manageable names
ovary.tib$Race <- ifelse(ovary.tib$Race == "Hispanic (All Races)", "Hisp",
        ifelse(ovary.tib$Race == "Non-Hispanic Asian or Pacific Islander", "API",
        ifelse(ovary.tib$Race == "Non-Hispanic American Indian/Alaska Native", "Native",
        ifelse(ovary.tib$Race == "Non-Hispanic Black", "Black",
        ifelse(ovary.tib$Race == "Non-Hispanic White", "White", "Unk")))))

#Change Nodes_Ex to None, Inadequate, Adequate
ovary.tib$Nodes_Ex <- ifelse(ovary.tib$Nodes_Ex == 0, "None",
        ifelse(ovary.tib$Nodes_Ex > 0 & ovary.tib$Nodes_Ex < 10, "Inadequate", "Adequate"))

#Change Nodes_Pos to Yes/No
ovary.tib$Nodes_Pos <- ifelse(ovary.tib$Nodes_Pos == 0, "No",
                    ifelse(ovary.tib$Nodes_Pos == 98, "Unk", "Yes"))

#Change COD variable to Alive(0) or Dead(1)
ovary.tib$COD <- ifelse(ovary.tib$COD == "Alive or dead of other cause", 0, 1)

#Filter out HGSOC and add column for black race and age groups
serous <- filter(ovary.tib, Histo == "8441" | Histo == "8460" | Histo == "8461" )
HGS <- filter(serous, Grade == 3 | Grade == 4)
HGS = mutate(HGS, Black.Race = ifelse(HGS$Race == "Black", "yes", "no") )
HGS$Age.Group <- cut(HGS$DiagAge, breaks = c(18, 30, 40, 50, 60 ,70, 80, 120), labels = c("18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"), right = FALSE )

#Filter those with HGSOC who got chemo
HGS.chemo <- filter(HGS, HGS$Chemo == "Yes")

#Filter HGSOC w/chemo and unknown nodal status (Nx)
HGS.chemoNx <- filter(HGS.chemo,HGS.chemo$Nodes_Ex == "None")

#Filter out those with negative nodes who received chemo all/24mo/36mo/60mo
HGS.chemo.N0removed <- filter(HGS.chemo, HGS.chemo$Nodes_Pos != "No")

HGS.chemo.N0removed24 <- filter(HGS.chemo.N0removed, HGS.chemo.N0removed$SurvMonths < 25)

HGS.chemo.N0removed36 <- filter(HGS.chemo.N0removed, HGS.chemo.N0removed$SurvMonths < 37)

HGS.chemo.N0removed60 <- filter(HGS.chemo.N0removed, HGS.chemo.N0removed$SurvMonths < 61)

HGS.chemo.N0removed120 <- filter(HGS.chemo.N0removed, HGS.chemo.N0removed$SurvMonths < 121)

#Filter HGSOC who did not get chemo/unknown chemo status
HGS.nochemo <- filter(HGS, HGS$Chemo != "Yes")

#Filter HGSOC who have Nx nodes overall/24mo/36mo/60mo
HGS.Nx <- filter(HGS, HGS$Nodes_Ex == "None")

HGS.Nx24 <- filter(HGS.Nx, HGS.Nx$SurvMonths < 25)

HGS.Nx36 <- filter(HGS.Nx, HGS.Nx$SurvMonths < 37)

HGS.Nx60 <- filter(HGS.Nx, HGS.Nx$SurvMonths < 61)

#Filter HGSOC who have positive nodes and received chemo
HGS.chemoN1 <- filter(HGS.chemo, HGS.chemo$Nodes_Pos == "Yes")


#Stratify HGSOC by stage and whether or not they had nodes evaluated
table(HGS$T_Stage, HGS$Nodes_Ex)
###    Adequate Inadequate None
###T1a      210        113  118
###T1b       81         37   22
###T1c      325        186  178

#Create KM survival fit line using library(survival) and library(survminer)


##fit <- survfit(Surv(time = [time variable], event = [censoring variable]) ~ [stratification variable], data = [dataset])
##ggsurvplot(fit, data = [dataset]) ##Lots of attributes available for custom plots##

#Fit equation for HGSOC stratified by LND all/with chemo/without chemo
fit.byLND <- survfit(Surv(time = HGS$SurvMonths, event = HGS$COD) ~ HGS$Nodes_Ex, data = HGS)

fit.byLND.chemo <- survfit(Surv(time = HGS.chemo$SurvMonths, event = HGS.chemo$COD) ~ HGS.chemo$Nodes_Ex, data = HGS.chemo)

fit.byLND.nochemo <- survfit(Surv(time = HGS.nochemo$SurvMonths, event = HGS.nochemo$COD) ~ HGS.nochemo$Nodes_Ex, data = HGS.nochemo)

#Fit equation for HGSOC receiving chemo stratified by nodal status
fit.byNodeStat.chemo <- survfit(Surv(time = HGS.chemo$SurvMonths, event = HGS.chemo$COD) ~ HGS.chemo$Nodes_Pos, data = HGS.chemo)

#Fit equation for HGSOC receiving chemo ~ nodal status with N0 removed all/24mo/36mo/60mo
fit.byNodeStatN0remove.chemo <- survfit(Surv(time = HGS.chemo.N0removed$SurvMonths, event = HGS.chemo.N0removed$COD) ~ HGS.chemo.N0removed$Nodes_Pos, data = HGS.chemo.N0removed)

fit.byNodeStatN0remove.chemo24 <- survfit(Surv(time = HGS.chemo.N0removed24$SurvMonths, event = HGS.chemo.N0removed24$COD) ~ HGS.chemo.N0removed24$Nodes_Pos, data = HGS.chemo.N0removed24)

fit.byNodeStatN0remove.chemo36 <- survfit(Surv(time = HGS.chemo.N0removed36$SurvMonths, event = HGS.chemo.N0removed36$COD) ~ HGS.chemo.N0removed36$Nodes_Pos, data = HGS.chemo.N0removed36)

fit.byNodeStatN0remove.chemo60 <- survfit(Surv(time = HGS.chemo.N0removed60$SurvMonths, event = HGS.chemo.N0removed60$COD) ~ HGS.chemo.N0removed60$Nodes_Pos, data = HGS.chemo.N0removed60)

fit.byNodeStatN0remove.chemo120 <- survfit(Surv(time = HGS.chemo.N0removed120$SurvMonths, event = HGS.chemo.N0removed120$COD) ~ HGS.chemo.N0removed120$Nodes_Pos, data = HGS.chemo.N0removed120)

#Fit equation for HGSOC receiving chemo ~ black race or not
fit.byRace <- survfit(Surv(time = HGS.chemo$SurvMonths, event = HGS.chemo$COD) ~ HGS.chemo$Black.Race, data = HGS.chemo)

#Survival curve for HGSOC stratified by LND with or without chemo
ggsurvplot(fit.byLND, data = HGS, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival stratified by LND with or without Chemo", legend = "bottom", legend.title = "LN status")

#Survival curve for HGSOC stratified by LND with chemo
ggsurvplot(fit.byLND.chemo, data = HGS.chemo, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival stratified by LND with Chemo", legend = "bottom", legend.title = "LN status")

#Survival curve for HGSOC stratified by LND without chemo
ggsurvplot(fit.byLND.nochemo, data = HGS.nochemo, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival stratified by LND with no/unknown Chemo", legend = "bottom", legend.title = "LN status")

#Survival curve for HGSOC receiving chemo ~ nodal status
ggsurvplot(fit.byNodeStat.chemo, data = HGS.chemo, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival in patients receiving Chemo ~ Nodal status", legend = "bottom", legend.title = "Nodal Status")

#Survival curve for HGSOC receiving chemo ~ nodal status with N0 removed all/24mo/36mo/60mo
ggsurvplot(fit.byNodeStatN0remove.chemo, data = HGS.chemo.N0removed, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival in patients receiving Chemo ~ Nodal status", legend = "bottom", legend.title = "Nodal Status")

ggsurvplot(fit.byNodeStatN0remove.chemo24, data = HGS.chemo.N0removed24, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival in patients receiving Chemo ~ Nodal status at 2 years", legend = "bottom", legend.title = "Nodal Status")

ggsurvplot(fit.byNodeStatN0remove.chemo36, data = HGS.chemo.N0removed36, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival in patients receiving Chemo ~ Nodal status at 3 years", legend = "bottom", legend.title = "Nodal Status")

ggsurvplot(fit.byNodeStatN0remove.chemo60, data = HGS.chemo.N0removed60, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival in patients receiving Chemo ~ Nodal status at 5 years", legend = "bottom", legend.title = "Nodal Status")

ggsurvplot(fit.byNodeStatN0remove.chemo120, data = HGS.chemo.N0removed120, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival in patients receiving Chemo ~ Nodal status at 10 years", legend = "bottom", legend.title = "Nodal Status")

#Survival curve for HGSOC receiving chemo ~ black race
ggsurvplot(fit.byRace.black, data = HGS.chemo.black, pval = TRUE)

#Fit equation for HGSOC with Nx nodes who received chemo overall
fit.Nx.Chemo <- survfit(Surv(time = HGS.chemoNx$SurvMonths, event = HGS.chemoNx$COD) ~ 1, data = HGS.chemoNx)

#Fit equation for HGSOC with N1 nodes who received chemo overall
fit.N1.Chemo <- survfit(Surv(time = HGS.chemoN1$SurvMonths, event = HGS.chemoN1$COD) ~ 1, data = HGS.chemoN1)

#Fit equation for HGSOC that did not receive chemo stratified by Nodal status/Race
fit.nochemo.byLND <- survfit(Surv(time = HGS.nochemo$SurvMonths, event = HGS.nochemo$COD) ~ HGS.nochemo$Nodes_Ex, data = HGS.nochemo)

fit.nochemo.byRace <- survfit(Surv(time = HGS.nochemo$SurvMonths, event = HGS.nochemo$COD) ~ HGS.nochemo$Black.Race, data = HGS.nochemo)
#Survival curve for HGSOC that did not receive chemo stratified by Nodal status/Race
ggsurvplot(fit.nochemo.byLND, data = HGS.nochemo, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival stratified by LND without Chemo", legend = "bottom", legend.title = "LN status")

ggsurvplot(fit.nochemo.byRace, data = HGS.nochemo, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival stratified by Race without Chemo", legend = "bottom", legend.title = "Race")

#Fit equation for HGSOC with Nx nodes stratified by chemo status overall/24mo/36mo/60mo
fit.Nx.byChemo <- survfit(Surv(time = HGS.Nx$SurvMonths, event = HGS.Nx$COD) ~ HGS.Nx$Chemo, data = HGS.Nx)

fit.Nx.byChemo24 <- survfit(Surv(time = HGS.Nx24$SurvMonths, event = HGS.Nx24$COD) ~ HGS.Nx24$Chemo, data = HGS.Nx24)

fit.Nx.byChemo36 <- survfit(Surv(time = HGS.Nx36$SurvMonths, event = HGS.Nx36$COD) ~ HGS.Nx36$Chemo, data = HGS.Nx36)

fit.Nx.byChemo60 <- survfit(Surv(time = HGS.Nx60$SurvMonths, event = HGS.Nx60$COD) ~ HGS.Nx60$Chemo, data = HGS.Nx60)

#Survplot for HGSOC with Nx nodes stratified by chemo status overall/36mo/60mo
ggsurvplot(fit.Nx.byChemo, data = HGS.Nx, pval = TRUE, xlab = "Months", break.time.by = 12, title = "Survival of Nx strat by Chemo", legend = "bottom", legend.title = "Chemo Received")

ggsurvplot(fit.Nx.byChemo24, data = HGS.Nx24, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival of Nx strat by Chemo at 2years", legend = "bottom", legend.title = "Chemo Received")

ggsurvplot(fit.Nx.byChemo36, data = HGS.Nx36, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival of Nx strat by Chemo at 3years", legend = "bottom", legend.title = "Chemo Received")

ggsurvplot(fit.Nx.byChemo60, data = HGS.Nx60, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival of Nx strat by Chemo at 5 years", legend = "bottom", legend.title = "Chemo Received")

##Combine more than one fit equation on survival plot

##Create list of fit equations
#fit.list <- list(Label1 = fit1, Label2 = fit2, etc.)

##Plot combined fit equation onto KM plot
#ggsurvplot_combine(fit.list, data = [dataset]) ##Lots of attributes available for custom plots##

#Create combined fit.list of HGSOC with Nx nodes that received chemo and HGSOC w/wo Chemo strat by LND
fit.combine.NXchemo.byLND <- list(Nx.Chemo = fit.Nx.Chemo, LND = fit.byLND)

#Create combined fit.list of HGSOC with N1 nodes that received chemo and HGSOC w/Chemo strat by LND
fit.combine.N1chemo.byLNDchemo <- list(N1 = fit.N1.Chemo, LND = fit.byLND.chemo)

#Plot HGSOC with Nx nodes that received chemo and HGSOC w/wo Chemo strat by LND
ggsurvplot_combine(fit.combine.NXchemo.byLND, data = HGS, pval = TRUE, xlab = "Months", break.time.by = 6, legend = "bottom")

#Plot HGSOC with N1 nodes on top of HGSOC stratified by LND with all patients receiving chemo
ggsurvplot_combine(fit.combine.N1chemo.byLNDchemo, data = HGS, risk.table = FALSE, ggtheme = theme_grey(), pval = TRUE, xlab = "Months", break.time.by = 6, legend = "bottom", title = "Positive nodes plotted on top of HGSOC ~ LND, all with chemo")


