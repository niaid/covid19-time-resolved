library(lme4)
library(lmerTest)

juncture.divergent.proteins <- function(data.matrix,pval.cutoff = 0.05) {
  selected.cytokine.profiles.permute <- data.matrix
  permuted.type1.proteins <- data.frame()
  permuted.type2.proteins <- data.frame()
  type1.model.count <- 0
  type2.model.count <- 0
  for (i in unique(selected.cytokine.profiles$test_id)) {
    if (sum(is.na(selected.cytokine.profiles.permute[,i])) < nrow(selected.cytokine.profiles.permute)/2) {

      tryCatch({
        # type 1
        # DE?
        # early stage reference
        form <- test_value ~ outcome + age + sex + days_from_symptom_onset_to_test + (1|subject_id)
        tmp <- subset(selected.cytokine.profiles.permute,stage == "early")
        tmp$test_value <- tmp[,i]
        lm.res <- lmer(form,tmp)
        lm.anova <- anova(lm.res)
        early.ref <- data.frame(coef=t(fixef(lm.res)),pval = t(lm.anova$`Pr(>F)`),singular=isSingular(lm.res))
        colnames(early.ref) <- paste0(colnames(early.ref),".early")
        
        for (s in c("mid","late")) {
          form <- test_value ~ outcome + age + sex + days_from_symptom_onset_to_test + (1|subject_id)
          tmp <- subset(selected.cytokine.profiles.permute,stage == s)
          tmp$test_value <- tmp[,i]#scale(log10(tmp[,i]+1))
          lm.res <- lmer(form,tmp)
          lm.anova <- anova(lm.res)
          de.res <- data.frame(coef=t(fixef(lm.res)),pval = t(lm.anova$`Pr(>F)`),singular=isSingular(lm.res),early.ref)
          type1.model.count <- type1.model.count + 1
          #if DE, test for significant change compared to early
          if (de.res$pval.1 <= pval.cutoff) {
            form <- test_value ~ test.interval*outcome + age + sex + days_from_symptom_onset_to_test + (1|subject_id)
            did.tmp <- subset(selected.cytokine.profiles.permute,stage %in% c("early",s))
            did.tmp$test.interval <- did.tmp$stage != "early"
            did.tmp$test_value <- did.tmp[,i]
            
            did.lm.res <- lmer(form,did.tmp)
            did.summary.lm <- summary(did.lm.res)
            
            did.df <- data.frame(coef=t(fixef(did.lm.res)),pval=t(did.summary.lm$coefficients[,"Pr(>|t|)"]),singular=isSingular(did.lm.res))
            if (did.df$pval.test.intervalTRUE.outcomeDeceased <= pval.cutoff)
              permuted.type1.proteins <- rbind(permuted.type1.proteins,data.frame(cytokine=i,stage=s,de=de.res,did=did.df))
          }
        }
        
        # type 2
        # DE between stages compared to early for each outcome
        # early.v.mid
        form <- test_value ~ test.interval + age + sex + (1|subject_id)
        type2.res <- data.frame()
        for (o in unique(selected.cytokine.profiles.permute$outcome)) {
          for (s in c("mid","late")) {
            tmp <- subset(selected.cytokine.profiles.permute,stage %in% c("early",s) & outcome == o)
            tmp$test.interval <- tmp$stage != "early"
            tmp$test_value <- tmp[,i]
            lm.res <- lmer(form,tmp)
            lm.anova <- anova(lm.res)
            test.res <- data.frame(coef=t(fixef(lm.res)),pval = t(lm.anova$`Pr(>F)`),singular=isSingular(lm.res))
            type2.res <- rbind(type2.res,data.frame(outcome=o,stage=s,test.res))
            type2.model.count <- type2.model.count + 1
          }
        }
        type2.res <- merge(subset(type2.res,outcome == "Deceased"),subset(type2.res,outcome == "Alive" & pval.1),
                           by="stage",suffixes=c(".dec",".rec"))
        type2.res <- subset(type2.res,pval.1.dec <= pval.cutoff & 
                              (pval.1.rec > pval.cutoff | sign(coef.test.intervalTRUE.rec) != sign(coef.test.intervalTRUE.dec)))
        if (nrow(type2.res) > 0)
          permuted.type2.proteins <- rbind(permuted.type2.proteins,data.frame(cytokine=i,type2.res))
        
      },error = function(e){})
    }
  }
  res <- list(type1=permuted.type1.proteins,type1.count=type1.model.count,type2=permuted.type2.proteins,type2.count=type2.model.count)
  return(res)
}