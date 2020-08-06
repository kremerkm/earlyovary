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
HGserous <- filter(serous, Grade == 3 | Grade == 4)

#Stratify HGSOC by stage and whether or not they had nodes evaluated
table(HGserous$T_Stage, HGserous$Nodes_Ex)
###    Adequate Inadequate None
###T1a      210        113  118
###T1b       81         37   22
###T1c      325        186  178

#Create KM survival fit line using library(survival) and library(survminer)
#
library(survival)

library(survminer)

##fit <- survfit(Surv(time = [time variable], event = [censoring variable]) ~ [stratification variable], data = [dataset])
#
###ggsurvplot(fitChemoNodes36, data = HGserousChemo36, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Chemo with nodal status at 3y", legend = "bottom", legend.title = "LN status")

#Create KM plot using survminer
#
##ggsurvplot(fit, data = [dataset]) ##Lots of attributes available for custom plots##
#
###ggsurvplot(fitChemoNodes36, data = HGserousChemo36, pval = TRUE, xlab = "Months", break.time.by = 6, title = "Chemo with nodal status at 3y", legend = "bottom", legend.title = "LN status")
