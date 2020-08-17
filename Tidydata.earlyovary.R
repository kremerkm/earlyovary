library(tidyverse)

library(SEER2R)
#Read SEER case files for ovarian cancer from 2004-2015 into data.frame
ovary.df <- read.SeerStat("ovcase_2004-2015.dic", UseVarLabelsInData = TRUE)

#Change data.frame to tibble
ovary.tib <- as_tibble(ovary.df)

#Rename columns
ovary.tib = rename(ovary.tib,
            Race = "Race_and_origin_recode_NHW_NHB_NHAIAN_NHAPI_Hispanic",
            Histo = "Histologic_Type_ICDO3",
            Lat = "Laterality",
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
            FollowMonth = "Month_of_followup_recode",
            FollowYear = "Year_of_followup_recode",
            CutoffStatus2 = "End_Calc_Vital_Status_Adjusted",
            MonthsFollowed = "Number_of_Intervals_Calculated")

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

#Combine Diagnosis Month and Diagnosis Year to Date and Follow-up Month/year to Date
ovary.tib = ovary.tib %>% 
        unite(DiagDate, c(DiagYear, DiagMonth), sep = "-")

library(zoo)
ovary.tib$DiagDate <- as.yearmon(ovary.tib$DiagDate, format = "%Y-%B")

ovary.tib = ovary.tib %>% 
        unite(FollowDate, c(FollowYear, FollowMonth), sep = "-")

ovary.tib$FollowDate <- as.yearmon(ovary.tib$FollowDate, format = "%Y-%B")

#Change COD variable to Alive(0) or Dead(1)
ovary.tib$COD <- ifelse(ovary.tib$COD == "Alive or dead of other cause", 0, 1)

#Filter out HGSOC
serous <- filter(ovary.tib, Histo == "441" | Histo == "460" | Histo == "461" )
HGS <- filter(serous, Grade == 3 | Grade == 4)

#Filter those with HGSOC who got chemo
HGS.chemo <- filter(HGS, HGS$Chemo == "Yes")

#Filter HGSOC w/chemo and unknown nodal status (Nx)
HGS.chemoNx <- filter(HGS.chemo,HGS.chemo$Nodes_Ex == "None")

#Filter out those with negative nodes who received chemo
HGS.chemo.N0removed <- filter(HGS.chemo, HGS.chemo$Nodes_Pos != "No")

#Filter HGSOC who did not get chemo/unknown chemo status
HGS.nochemo <- filter(HGS, HGS$Chemo != "Yes")

#Filter HGSOC who have Nx nodes overall/36mo/60mo
HGS.Nx <- filter(HGS, HGS$Nodes_Ex == "None")

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

library(survival)

library(survminer)

##fit <- survfit(Surv(time = [time variable], event = [censoring variable]) ~ [stratification variable], data = [dataset])
##ggsurvplot(fit, data = [dataset]) ##Lots of attributes available for custom plots##

#Fit equation for HGSOC stratified by LND all/with chemo/without chemo
fit.byLND <- survfit(Surv(time = HGS$SurvMonths, event = HGS$COD) ~ HGS$Nodes_Ex, data = HGS)

fit.byLND.chemo <- survfit(Surv(time = HGS.chemo$SurvMonths, event = HGS.chemo$COD) ~ HGS.chemo$Nodes_Ex, data = HGS.chemo)

fit.byLND.nochemo <- survfit(Surv(time = HGS.nochemo$SurvMonths, event = HGS.nochemo$COD) ~ HGS.nochemo$Nodes_Ex, data = HGS.nochemo)

#Fit equation for HGSOC receiving chemo stratified by nodal status
fit.byNodeStat.chemo <- survfit(Surv(time = HGS.chemo$SurvMonths, event = HGS.chemo$COD) ~ HGS.chemo$Nodes_Pos, data = HGS.chemo)

#Fit equation for HGSOC receiving chemo ~ nodal status with N0 removed
fit.byNodeStatN0remove.chemo <- survfit(Surv(time = HGS.chemo.N0removed$SurvMonths, event = HGS.chemo.N0removed$COD) ~ HGS.chemo.N0removed$Nodes_Pos, data = HGS.chemo.N0removed)

#Survival curve for HGSOC stratified by LND with or without chemo
ggsurvplot(fit.byLND, data = HGS, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival stratified by LND with or without Chemo", legend = "bottom", legend.title = "LN status")

#Survival curve for HGSOC stratified by LND with chemo
ggsurvplot(fit.byLND.chemo, data = HGS.chemo, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival stratified by LND with Chemo", legend = "bottom", legend.title = "LN status")

#Survival curve for HGSOC stratified by LND without chemo
ggsurvplot(fit.byLND.nochemo, data = HGS.nochemo, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival stratified by LND with no/unknown Chemo", legend = "bottom", legend.title = "LN status")

#Survival curve for HGSOC receiving chemo ~ nodal status
ggsurvplot(fit.byNodeStat.chemo, data = HGS.chemo, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival in patients receiving Chemo ~ Nodal status", legend = "bottom", legend.title = "Nodal Status")

#Survival curve for HGSOC receiving chemo ~ nodal status with N0 removed
ggsurvplot(fit.byNodeStatN0remove.chemo, data = HGS.chemo.N0removed, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Survival in patients receiving Chemo ~ Nodal status", legend = "bottom", legend.title = "Nodal Status")

#Fit equation for HGSOC with Nx nodes who received chemo overall
fit.Nx.Chemo <- survfit(Surv(time = HGS.chemoNx$SurvMonths, event = HGS.chemoNx$COD) ~ 1, data = HGS.chemoNx)

#Fit equation for HGSOC with N1 nodes who received chemo overall
fit.N1.Chemo <- survfit(Surv(time = HGS.chemoN1$SurvMonths, event = HGS.chemoN1$COD) ~ 1, data = HGS.chemoN1)

#Fit equation for HGSOC that did not receive chemo stratified by Nodal status
fit.nochemo.byLND <- survfit(Surv(time = HGS.nochemo$SurvMonths, event = HGS.nochemo$COD) ~ HGS.nochemo$Nodes_Ex, data = HGS.nochemo)

#Survival curve for HGSOC that did not receive chemo stratified by Nodal status at 5 years
ggsurvplot(fit.nochemo.byLND, data = HGS.nochemo, pval = TRUE, xlab = "Months", break.time.by = 6, xlim = c(0,60), title = "Survival stratified by LND without Chemo", legend = "bottom", legend.title = "LN status")

#Fit equation for HGSOC with Nx nodes stratified by chemo status overall/36mo/60mo
fit.Nx.byChemo <- survfit(Surv(time = HGS.Nx$SurvMonths, event = HGS.Nx$COD) ~ HGS.Nx$Chemo, data = HGS.Nx)

fit.Nx.byChemo36 <- survfit(Surv(time = HGS.Nx36$SurvMonths, event = HGS.Nx36$COD) ~ HGS.Nx36$Chemo, data = HGS.Nx36)

fit.Nx.byChemo60 <- survfit(Surv(time = HGS.Nx60$SurvMonths, event = HGS.Nx60$COD) ~ HGS.Nx60$Chemo, data = HGS.Nx60)

#Survplot for HGSOC with Nx nodes stratified by chemo status overall/36mo/60mo
ggsurvplot(fit.Nx.byChemo, data = HGS.Nx, pval = TRUE, xlab = "Months", break.time.by = 12, title = "Survival of Nx strat by Chemo", legend = "bottom", legend.title = "Chemo Received")

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


