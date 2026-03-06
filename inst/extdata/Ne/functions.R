#================================================================================================
# 格式: utf-8
# 版权: GPLv3 英文 https://www.gnu.org/licenses/quick-guide-gplv3.en.html
#             中文 https://www.gnu.org/licenses/quick-guide-gplv3.html
# 作者：栾生, luansheng@ysfri.ac.cn
# 功能: 定义一些有用的小函数，主要包括：
#       1.读取ASReml单机版 4.1 输出的xml文件，获得育种值、方差组分、均值等参数；
#       2.预测和实现遗传增益的计算函数；
#       3.表型数据结构分析，包括全半同胞家系数量和比例、父母本数量等；
#       4.描述性统计参数计算函数。
#================================================================================================

# 从ASReml 4.1的xml输出文件中，获得育种值。
# 使用R包xml2解析xml文件。
readEBV <- function(asreml.xml, EBV.term) {
  #\\s+正则表达式，其中第一个\表示开始正则表达式，
  #\s 表示空格，+表示1个或者多个空格
  # solutions <-
  #   as.data.table(tstrsplit(trimws(xml_text(
  #     xml_find_all(asreml.xml, "/ASReport/Cycle/Solutions/Equation")
  #   )), "\\s+"))
  ModelTerm <- xml_text(
    xml_find_all(asreml.xml, "/ASReport/Cycle/Solutions/Equation/Model_Term")
  )
  Level <- xml_text(
    xml_find_all(asreml.xml, "/ASReport/Cycle/Solutions/Equation/Level")
  )
  Effect <- xml_text(
    xml_find_all(asreml.xml, "/ASReport/Cycle/Solutions/Equation/Effect")
  )
  seEffect <- xml_text(
    xml_find_all(asreml.xml, "/ASReport/Cycle/Solutions/Equation/seEffect")
  )
  solutions <- data.table(ModelTerm,Level,Effect,seEffect)
  setnames(solutions,
           colnames(solutions),
           c("ModelTerm", "TraitID", "Effect", "EffectSE"))

  # 基于A矩阵，提取单性状父母本模型的解
  if (EBV.term == "Sire.Dam") {
    sire.dam.effects <-
      solutions[ModelTerm == "nrmv(SireID)", 1:4][, 2:4]
    colnames(sire.dam.effects) <-
      c("AnimalID", "Effect", "EffectSE")
    sire.dam.effects[, ":="(Effect = as.numeric(Effect),
                            EffectSE = as.numeric(EffectSE))]
    return(sire.dam.effects)
  }
  
  # 基于H矩阵，提取单性状父母本模型的解
  if (EBV.term == "Sire.DamH") {
    sire.dam.effects <-
      solutions[ModelTerm == "giv1(SireID)", 1:4][, 2:4]
    colnames(sire.dam.effects) <-
      c("AnimalID", "Effect", "EffectSE")
    sire.dam.effects[, ":="(Effect = as.numeric(Effect),
                            EffectSE = as.numeric(EffectSE))]
    return(sire.dam.effects)
  }
  # #ASReml4 新格式 读取cycle格式xml数据，貌似有问题
  # if (EBV.term == "Cycle.Animal") {
  #   animal.effects <-
  #     solutions[ModelTerm == "nrmv(AnimalID)", 1:4][, 2:4]
  #   colnames(animal.effects) <- c("AnimalID", "Effect", "EffectSE")
  #   animal.effects[, Effect := as.numeric(Effect)]
  #   animal.effects[, EffectSE := as.numeric(EffectSE)]
  #   return(animal.effects)
  # }
  # 基于A矩阵，提取单性状个体动物模型的解
  if (EBV.term == "Animal") {
    animal.effects <-
      solutions[ModelTerm == "nrmv(AnimalID)", 1:4][, 2:4]
    colnames(animal.effects) <- c("AnimalID", "Effect", "EffectSE")
    animal.effects[, Effect := as.numeric(Effect)]
    animal.effects[, EffectSE := as.numeric(EffectSE)]
    return(animal.effects)
  }

  # 基于H矩阵，提取单性状个体动物模型的解
  if (EBV.term == "AnimalH") {
    animal.effects <-
      solutions[ModelTerm == "giv1(AnimalID)", 1:4][, 2:4]
    colnames(animal.effects) <- c("AnimalID", "Effect", "EffectSE")
    animal.effects[, Effect := as.numeric(Effect)]
    animal.effects[, EffectSE := as.numeric(EffectSE)]
    return(animal.effects)
  }

  # 基于A矩阵，提取两性状个体动物模型的解
  if (EBV.term == "Tr.Animal") {
    animal.effects <-
      solutions[ModelTerm == "us(Tr).nrm(AnimalID)", 1:4][, 2:4]
    animal.effects[, ":="(TraitID = substr(TraitID, 1, 1),
                          AnimalID = substring(TraitID, 3))]
    animal.effects[, ":="(
      TraitID = as.numeric(TraitID),
      Effect = as.numeric(Effect),
      EffectSE = as.numeric(EffectSE)
    )]
    return(animal.effects)
  }
  
  # 基于H矩阵，提取两性状个体动物模型的解
  if (EBV.term == "Tr.AnimalH") {
    animal.effects <-
      solutions[ModelTerm == "us(Tr).giv1(AnimalID)", 1:4][, 2:4]
    animal.effects[, ":="(TraitID = substr(TraitID, 1, 1),
                          AnimalID = substring(TraitID, 3))]
    animal.effects[, ":="(
      TraitID = as.numeric(TraitID),
      Effect = as.numeric(Effect),
      EffectSE = as.numeric(EffectSE)
    )]
    return(animal.effects)
  }

  # 提取总体均值
  if (EBV.term == "mu") {
    animal.effects <- solutions[ModelTerm == "mu", 1:4][, 2:4]
    colnames(animal.effects) <- c("AnimalID", "Effect", "EffectSE")
    animal.effects[, Effect := as.numeric(Effect)]
    animal.effects[, EffectSE := as.numeric(EffectSE)]
    return(animal.effects)
  }
}

# 定义一个模板函数，用于生成提取育种值的多个函数
read_env_ebv <- function(model_label) {
  function(env_seq, env_asreml_xml_filename) {
    tmp_list <- vector("list", length = nrow(env_seq))
    for (i in 1:nrow(env_seq)) {
      env_asreml_xml <-
        read_xml(env_asreml_xml_filename[i], encoding = "UTF-8")
      env_effects <-
        readEBV(asreml.xml = env_asreml_xml, EBV.term = model_label)
      env_effects[, ":="(EnvID = rep(env_seq$EnvID[i], .N))]
      tmp_list[[i]] <- env_effects
    }
    return(rbindlist(tmp_list, use.names = TRUE))
  }
}

# 基于模板函数，定义新的函数，读取不同模型的解
# 调用新函数时，需要读入env_seq和env_asreml_xml_filename两个参数
# 基于A矩阵，读取单性状个体动物模型的EBV
read_env_ebv_animal <- read_env_ebv("Animal")

# 基于H矩阵，读取单性状个体动物模型的EBV
read_env_ebv_animal_H <- read_env_ebv("AnimalH")

# 基于A矩阵，读取单性状父母本模型的EBV
read_env_ebv_siredam <- read_env_ebv("Sire.Dam")

# 基于H矩阵，读取单性状父母本模型的EBV
read_env_ebv_siredam_H <- read_env_ebv("Sire.DamH")

## 基于H矩阵，读取两性状个体动物模型的EBV
read_env_ebv_trait_animal_H <- read_env_ebv("Tr.AnimalH")

# 读取单性状模型的群体均值
read_env_ebv_mu <- read_env_ebv("mu")



# # 从ASReml 4.1的xml输出文件中，获得方差组分。
# 使用R包xml2解析xml文件。
readvc <- function(asreml.xml) {
  variance.components.source <-
    trimws(xml_text(
      xml_find_all(
        asreml.xml,
        "/ASReport/Cycle/VarianceComponents/VParameter/Source"
      )
    ))
  variance.components.value <-
    as.numeric(xml_text(
      xml_find_all(
        asreml.xml,
        "/ASReport/Cycle/VarianceComponents/VParameter/VComponent"
      )
    ))
  g.vc.values <-
    data.table(variance.components.source, variance.components.value)
  return(g.vc.values)
}

# 定义模板函数，读取不同模型的方差组分
# 如果是多性状模型，譬如两性状模型，会输出多个性状的方差和协方差
read_env_avc <- function(model_label) {
  function(env_seq, env_asreml_xml_filename) {
    tmp_list <- vector("list", length = nrow(env_seq))
    for (i in 1:nrow(env_seq)) {
      env_asreml_xml <-
        read_xml(env_asreml_xml_filename[i], encoding = "UTF-8")
      env_vc_value <- readvc(env_asreml_xml)
      env_vc_value[, ":="(EnvID = rep(env_seq$EnvID[i], .N))]
      tmp_list[[i]] <- env_vc_value
    }
    vc_value <- rbindlist(tmp_list, use.names = TRUE)
    vc_value_additive <-
      vc_value[grepl(model_label, variance.components.source, fixed =
                       TRUE), variance.components.value]
    times <- length(vc_value_additive) / nrow(env_seq)
    # env_seq[rep(seq_len(nrow(env_seq)), each=times)]重复数据框中的每一行times次
    env_seq_avc <- data.table(env_seq[rep(seq_len(nrow(env_seq)), each=times)], aVC = vc_value_additive)
    return(env_seq_avc)
  }
}

# 基于模板函数，定义新的函数，读取不同模型的方法组分
# 调用新函数时，需要读入env_seq和env_asreml_xml_filename两个参数
# 基于A矩阵，读取单性状个体动物模型的加性方差
read_env_gavc <- read_env_avc(c("nrmv(AnimalID)"))
# 基于H矩阵，读取单性状个体动物模型的加性方差
read_env_gavc_H <- read_env_avc(c("giv1(AnimalID)"))

# 基于H矩阵，读取两性状个体动物模型的加性方差
read_env_gavc_trait_H <-  read_env_avc(c("us(Tr).giv1(AnimalID)"))

# 基于A矩阵，读取单性状父母本模型的加性方差
read_env_savc <- read_env_avc(c("nrmv(SireID)"))

# 基于H矩阵，读取单性状父母本模型的加性方差
read_env_savc_H <- read_env_avc(c("giv1(SireID)"))



# 提取近交系数的函数
read_env_f <- function(env_seq, filename) {
  tmp_list <- vector("list", length = nrow(env_seq))
  for (i in 1:nrow(env_seq)) {
    f <-
      fread(filename[i],
            header = T,
            stringsAsFactors = F)
    f[, ":="(EnvID = rep(env_seq$EnvID[i], .N))]
    tmp_list[[i]] <- f
  }
  return(rbindlist(tmp_list, use.names = TRUE))
}


# 预测和实现遗传进展计算
# 数据集中最少要包括2代的数据，必须要包括的字段：世代Generation、养殖环境EnvLabel、家系编号FamilyID。
# 计算思路：
# 不同世代收获体重的育种值差值，即为预测遗传进展
# 以体重为例，首先提取要计算世代的个体数据，然后获得他们的育种值。然后计算当前世代的育种值均值，据此在下一世代家系中寻找育种值与此最为相近的一定数量的家系，设置作为对照群体。
# 具体操作，思路是这样的，把要设置对找家系的世代所有家系的育种值，从小到大排列，然后按照递增顺序，按照对照家系数量一组，计算每组的均值，然后与当前世代的均值进行比较，差异最小的组，作为对照家系。
calc_sel_rp <-
  function(data,
           #数据集
           ct_fam_num,
           #对照家系数量
           tr_name,
           #待分析的性状字段名称
           tr_ebv_name,
           #性状育种值字段名称
           lmer_model) {
    #计算最小二乘均值的lmer模型
    ct_seq <- seq(ct_fam_num)
    gen_name <- unique(data$Generation)
    gen_num <- length(gen_name)
    env_label <- unique(data$EnvLabel)
    env_num <- length(env_label)
    if (("Generation" %in% colnames(data)) &
        (gen_num > 1)) {
      #计算预测选择反应
      gen_env_ebv <- data[, lapply(.SD, mean, na.rm = TRUE),
                          by = .(Generation, EnvLabel),
                          .SDcols = c(tr_name,
                                      tr_ebv_name)]
      gen_fam_num_size <- g_s_data_ebv_inbred[,
                                              .(.N,
                                                FamilyNum = length(unique(FamilyID))),
                                              by = .(Generation, EnvLabel)]
      gen_env_ebv <- merge(
        gen_env_ebv,
        gen_fam_num_size,
        by = c("Generation", "EnvLabel"),
        all.x = TRUE
      )
      #计算实现选择反应
      gen_fam_ebv <-
        data[, lapply(.SD, mean, na.rm = TRUE),
             by = .(Generation, EnvID, EnvLabel, FamilyID),
             .SDcols = c(tr_name, tr_ebv_name)]
      #排序该函数可以接受外部自定义字符串,全部为升序
      setorderv(gen_fam_ebv,
                c("Generation", "EnvLabel", tr_ebv_name),
                c(1, 1, 1))
      #给数据框加一个标识列
      gen_fam_ebv[, ":="(GenerationEnv = paste(Generation, EnvLabel, sep = ""))]
      #给数据框加一个标识列
      gen_env_ebv[, ":="(GenerationEnv = paste(Generation, EnvLabel,  sep = ""))]
      #把所有家系类型,首先设置为SP
      gen_fam_ebv[, ":="(SPCP = rep("SP", .N))]
      
      for (i in 1:(gen_num - 1)) {
        for (j in 1:env_num) {
          gen_env_name <-
            paste(gen_name[i], env_label[j], sep = "")
          next_gen_env_name <-
            paste(gen_name[i + 1], env_label[j], sep = "")
          #获得当前世代性状育种值均值
          cu_gen_ebv <-
            gen_env_ebv$gAnimalEBV[gen_env_ebv$GenerationEnv == gen_env_name]
          #获得下一世代家系育种值
          next_gen_fam_ebv <-
            gen_fam_ebv[(GenerationEnv %chin% next_gen_env_name)][order(Generation,
                                                                        EnvLabel,
                                                                        gAnimalEBV)]
          fam_num <- nrow(next_gen_fam_ebv)
          #如果下一代有100个家系，对照家系15个，那么就会有85个组，每组包括连续的15个家系
          #譬如1:15, 2:16, 3:17, 4:18等
          each_ct_mean <- array(0, fam_num - ct_fam_num)
          each_ct_seq <- ct_seq
          for (k in 1:(fam_num - ct_fam_num)) {
            each_ct_mean[k] <- mean(next_gen_fam_ebv[each_ct_seq, gAnimalEBV])
            #生成下一组，譬如从1:15 到 2:16
            each_ct_seq <- ct_seq + k
          }
          
          #寻找与当前世代育种值均值cu_gen_ebv最为接近的一个组
          #并获取家系在数据集中的位置，给这些家系打上CP标记
          ct_fam_pos <-
            ct_seq + which.min(abs(each_ct_mean - cu_gen_ebv))
          ct_family <- next_gen_fam_ebv[ct_fam_pos, FamilyID]
          gen_fam_ebv[(GenerationEnv %chin% next_gen_env_name) &
                        (FamilyID %chin% ct_family), ":="(SPCP = rep("CP", .N))]
        }
      }
    }
    #调用lmmer和emmeans两个包计算最小二乘均值
    #合并数据文件与育种值文件
    data_ebv_spcp <-
      merge(
        data,
        gen_fam_ebv[, .(Generation, EnvLabel, FamilyID, SPCP)],
        by = c("Generation", "EnvLabel", "FamilyID"),
        all.x = TRUE
      )
    data_ebv_spcp[, ":="(GenerationEnv = paste(Generation, EnvLabel, sep = ""))]
    
    gen_env_lmer_emmeans <-
      vector("list", (gen_num - 1) * env_num)
    
    #SP和CP两个值需要两个label，所以需要乘以2
    gen_env_label <-
      vector("character", 2 * (gen_num - 1) * env_num)
    k <- 0
    l <- 1
    for (i in 2:gen_num) {
      for (j in 1:env_num) {
        k <- k + 1
        gen_env_name <-
          paste(gen_name[i], env_label[j], sep = "")
        gen_env_label[l:(l + 1)] <-
          rep(gen_env_name, 2)
        gen_env_subset <-
          data_ebv_spcp[GenerationEnv %chin% gen_env_name,]
        gen_env_lmer <-
          eval(parse(
            text = paste("lmer(",
                         lmer_model,
                         ", data = gen_env_subset)",
                         sep = "")
          ))
        gen_env_lmer_emmeans[[k]] <-
          summary(emmeans(gen_env_lmer, ~ SPCP))
        #这里的2对应的是需要SP和CP群体2个label
        l <- l + 2
      }
    }
    gen_env_lmer_emmeans_dt <-
      as.data.table(cbind(GenerationEnv = gen_env_label,
                          do.call("rbind", gen_env_lmer_emmeans)))
    fam_num_per_gen_spcp <-
      data_ebv_spcp[, .(.N, FamilyNum = length(unique(FamilyID))), by = .(GenerationEnv, SPCP)]
    lsmean_per_gen_spcp <-
      merge(
        gen_env_lmer_emmeans_dt,
        fam_num_per_gen_spcp,
        by = c("GenerationEnv", "SPCP"),
        all.x = TRUE
      )
    return(
      list(
        PredictedSR = gen_env_ebv[order(EnvLabel, Generation)],
        RealizedSR = lsmean_per_gen_spcp[order(GenerationEnv, SPCP)],
        SPCPData = data_ebv_spcp
      )
    )
}

# 计算个体的选择指数
scale_100 <- function(x) {
  return(100 * (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}
# 表型数据结构分析，包括个体数量、全（半）同胞家系数量和比例、父母本数量等
summaryfamilystructure <-  function(individual.data) {
  individual.data.table <- as.data.table(individual.data)
  #个体的数量
  animal.num <- individual.data.table[, length(unique(AnimalID))]
  #亲本的数量
  sire.num <-   individual.data.table[, length(unique(SireID))]
  dam.num <- individual.data.table[, length(unique(DamID))]
  #分组统计每个家系个体的数量
  family.animal.num <-
    individual.data.table[, .N, by = .(SireID, DamID, FamilyID)]
  #分组统计与每个父本交配的母本数
  dam.num.per.sire <- family.animal.num[, .N, by = .(SireID)]
  
  #animal.num.per.fullsibfamily <- family.animal.num[, .(FamilyID, N)]
  #setnames(animal.num.per.fullsibfamily, c("N"), c("animal.num"))
  
  #全同胞家系的数量
  fullsib.family.num <- family.animal.num[, .N]
  
  #获得父系半同胞家系组的数量
  paternal.halfsibgroup.num <- dam.num.per.sire[N > 1, .N]
  #获得父系半同胞家系的数量
  paternal.halfsibfamily.num <- dam.num.per.sire[N > 1, sum(N)]
  
  #标记父系半同胞家系
  family.animal.num[, is.halfsib.family := array(FALSE, nrow(family.animal.num))]
  family.animal.num[, is.halfsib.family := is.halfsib.family |
                      (SireID %in% as.character(dam.num.per.sire[N > 1, SireID]))]
  
  #  sire.num.per.dam <-
  #    aggregate(cbind(SireID) ~ DamID, data = family.animal.num, FUN = length)
  sire.num.per.dam <- family.animal.num[, .N, by = .(DamID)]
  #获得母系半同胞家系组的数量
  maternal.halfsibgroup.Num <- sire.num.per.dam[N > 1, .N]
  #获得母系半同胞家系的数量
  maternal.halfsibfamily.num <-
    sire.num.per.dam[N > 1, sum(N)]
  #标记母系半同胞家系
  family.animal.num[, is.halfsib.family := is.halfsib.family |
                      (DamID %in% as.character(sire.num.per.dam[N > 1, DamID]))]
  
  #父系和母系半同胞家系合并的数量
  halfsibfamily.num <-
    family.animal.num[is.halfsib.family == TRUE, .N]
  #
  family.summary <-
    data.frame(
      AnimalNum = animal.num,
      SireNum = sire.num,
      DamNum = dam.num,
      FullsibFamilyNum = fullsib.family.num,
      HalfsibFamilyNum = halfsibfamily.num,
      PaternalHalfsibGroupNum = paternal.halfsibgroup.num,
      PaternalHalfsibFamilyNum = paternal.halfsibfamily.num,
      MaternalHalfsibGroupNum = maternal.halfsibgroup.Num,
      MaternalHalfsibFamilyNum = maternal.halfsibfamily.num
    )
  return(family.summary)
}

# 描述性统计函数
desc_stats <- function(x) {
  funs <- c(mean, min, max, sd, cv)
  if (sum(is.na(x)) == length(x)) {
    return(as.double(rep(NA, length(funs))))
  } else {
    return(as.double(unlist(lapply(funs, function(f)
      f(x, na.rm = TRUE)))))
  }
}

# 计算变异系数
cv <- function(x, na.rm = TRUE) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}


AnalyzeFamilyStructure <-  function(FamilyData) {
  FamilyData <- unique(FamilyData[,c("FamilyID","SireID","DamID")])
  FamilyData_tbl <- tbl_df(FamilyData[,c("FamilyID","SireID","DamID")])
  #全同胞家系的数量
  TotalFullSibFamilyNum <- nrow(FamilyData_tbl)
  #统计与每个父本交配的母本数目
  DamNumperSire <- tally(group_by(FamilyData_tbl, SireID))
  SireNum <- nrow(DamNumperSire)
  #获得父系半同胞家系组
  PaternalHalfSibGroup <- filter(DamNumperSire, n > 1)
  #获得父系半同胞家系组的数量
  PaternalHalfSibGroupNum <- nrow(PaternalHalfSibGroup)
  #获得父系半同胞家系的数量
  PaternalHalfSibFamilyNum <- sum(PaternalHalfSibGroup$n)
  
  #标记父系半同胞家系
  IsHalfSibFamily <- array(FALSE, TotalFullSibFamilyNum)
  FamilyData <- mutate(FamilyData, IsHalfSibFamily)
  FamilyData$IsHalfSibFamily <-
    FamilyData$IsHalfSibFamily |
    (FamilyData$SireID %in% as.character(PaternalHalfSibGroup$SireID))
  
  #统计与每个母本交配的父本数目
  SireNumperDam <- tally(group_by(FamilyData_tbl, DamID))
  DamNum <- nrow(SireNumperDam)
  #获得母系半同胞家系组
  MaternalHalfSibGroup <- filter(SireNumperDam, n > 1)
  #获得母系半同胞家系组的数量
  MaternalHalfSibGroupNum <- nrow(MaternalHalfSibGroup)
  #获得母系半同胞家系的数量
  MaternalHalfSibFamilyNum <- sum(MaternalHalfSibGroup$n)
  #标记母系半同胞家系
  FamilyData$IsHalfSibFamily <-
    FamilyData$IsHalfSibFamily |
    (FamilyData$DamID %in% as.character(MaternalHalfSibGroup$DamID))
  
  #父系和母系半同胞家系合并的数量
  TotalHalfSibFamilyNum <-
    length(FamilyData$FamilyID[FamilyData$IsHalfSibFamily])
  FamilySummary <-
    list(
      SireNum = SireNum,
      DamNum = DamNum,
      TotalFullSibFamilyNum = TotalFullSibFamilyNum,
      PaternalHalfSibGroupNum = PaternalHalfSibGroupNum,
      PaternalHalfSibFamilyNum = PaternalHalfSibFamilyNum,
      MaternalHalfSibGroupNum = MaternalHalfSibGroupNum,
      MaternalHalfSibFamilyNum = MaternalHalfSibFamilyNum,
      TotalHalfSibFamilyNum = TotalHalfSibFamilyNum,
      HalfsibRatio = sprintf("%.2f", TotalHalfSibFamilyNum * 100 / TotalFullSibFamilyNum)
    )
  return(FamilySummary)
}

# 生成类似facet_wrap( ~ Generation + SexID) 形式
ParseFacetTerms <- function(nested.terms = "",scales.value="fixed"){
  return(eval(parse(text=paste("facet_wrap(~",paste(unlist(strsplit(nested.terms,split = ",")),collapse = "+"),",scales=\"",scales.value,"\"",")",sep=""))))
}

#家系生产同步性分布柱状图
family_birthday_bar_theme <- theme_gray(base_size = 18,base_family = "heiti")+theme(
  axis.title.x=element_text(size=15),
  axis.title.y=element_text(size=15),
  axis.text.x=element_text(size=13,angle = 45, hjust = 1),
  axis.text.y=element_text(size=13),
  axis.ticks = element_line(size = 0.1),
  legend.position = "none",
  panel.border = element_blank()
)

SummarizeTableTrait <-
  function(trait.data,
           terms,
           terms.names,
           trait.column.name,
           trait.unit = "",
           count.unit = "",
           digit.num = 0,
           title) {
    if (trait.unit == "%") {
      eval(parse(text=paste("trait.data$",trait.column.name," <- trait.data$",trait.column.name," * 100",sep="")))
    }
    trait.data.groups <-
      eval(parse(text = paste(
        "group_by(trait.data,", terms, ")", sep = ""
      ))) 
    trait.data.groups.summ <-
      eval(parse(
        text = paste(
          "summarise(trait.data.groups,T.count=length(",
          trait.column.name,
          "),T.mean = mean(",
          trait.column.name,
          ",na.rm=TRUE),T.min = min(",
          trait.column.name,
          ",na.rm=TRUE),T.max=max(",
          trait.column.name,
          ",na.rm=TRUE),T.sd=sd(",
          trait.column.name,
          ",na.rm=TRUE),T.cv=cv(",
          trait.column.name,
          ")*100)",
          sep = ""
        )
      ))
    
    if (nrow(trait.data.groups.summ) > 1) {
      #汇总数据总体
      trait.data.summ <-
        eval(parse(
          text = paste(
            "summarise(trait.data,T.count=length(",
            trait.column.name,
            "),T.mean = mean(",
            trait.column.name,
            ",na.rm=TRUE),T.min=min(",
            trait.column.name,
            ",na.rm=TRUE),T.max=max(",
            trait.column.name,
            ",na.rm=TRUE),T.sd=sd(",
            trait.column.name,
            ",na.rm=TRUE),T.cv=cv(",
            trait.column.name,
            ")*100)",
            sep = ""
          )
        ))
      trait.data.groups.all.summ <-
        bind_rows(trait.data.groups.summ, trait.data.summ)
    } else {
      trait.data.groups.all.summ <- trait.data.groups.summ
    }
    
    Trait.kable <-
      kable(
        trait.data.groups.all.summ,
        col.names = c(unlist(strsplit(terms.names, ",")), paste("数量/",count.unit,sep=""), paste(c("平均值/", "最小值/", "最大值/", "标准差/"),trait.unit,sep=""), "变异系数/%"),
        digits = digit.num,
        caption = title
      )
    return(list(
      data = trait.data,
      summ = data.frame(trait.data.groups.all.summ),
      table = Trait.kable
    ))
  }

DescStatChinese <- function(vector.desc.stat,scale.unit="",digit.num=2){
  return(paste(
    "平均值为",
    sprintf(paste("%.",digit.num,"f",sep=""),mean(vector.desc.stat,na.rm=TRUE)),scale.unit,
    "，",
    "最小值为",
    sprintf(paste("%.",digit.num,"f",sep=""),min(vector.desc.stat,na.rm=TRUE)),scale.unit,
    "，",
    "最大值为",
    sprintf(paste("%.",digit.num,"f",sep=""),max(vector.desc.stat,na.rm=TRUE)),scale.unit,
    "，",
    "标准差为",
    sprintf(paste("%.",digit.num,"f",sep=""),sd(vector.desc.stat,na.rm=TRUE)),scale.unit,
    "，",
    "变异系数为",
    sprintf(paste("%.",digit.num,"f",sep=""),cv(vector.desc.stat)*100),
    " %",
    sep = ""
  )
  )
  
}  