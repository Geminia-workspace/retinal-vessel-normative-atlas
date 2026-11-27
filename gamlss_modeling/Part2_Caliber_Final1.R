# ========================================================================================
# Part 3: Caliber（口径，41 个参数）GAMLSS 分析  v15.1
#
# 【v15.1 核心更新】
#   1. 新 6 列分段布局：41 个参数分为 4 个小节，每节带中文小标题
#      - 动脉各分级平均管径（10）/ 静脉各分级平均管径（10）
#      - 1-3 级前 1/3/6 根平均直径（18）/ AVR（3）
#   2. 参数名简洁化：
#      artery1_avg_diameter_um  → "Diameter Artery_L1"
#      artery1_diameter_top1_um → "Diameter Artery_L1_Top1"
#      AVR_1                    → "AVR_L1"
#   3. Fig4a_Selected_Summary_MultiDisease 重写为分段布局（与 Trajectories 大图一致）
#      - 默认包含全部 41 参数；每节内部可通过 SECTIONS[[i]]$selected 限定子集
#   4. Fig5 新设计：显著性热图（替代旧的相关性热图）
#      - 行 = 参数，列 = 疾病 × 性别（CKD-F / CKD-M / AD-F / AD-M）
#      - 颜色 = median(Percentile) - 50（红=偏高, 蓝=偏低）
#      - 文本 = 中位数 + 显著性星号（Wilcoxon vs P50）
#   5. file_tag 不变 → 现有 RDS 无需重建；加载时自动用新 label 覆盖
#   6. 保留 ensure_avr()：当数据缺少 AVR_1/2/3 列时自动从 top6 直径计算
# ========================================================================================

suppressPackageStartupMessages({
  library(gamlss,     quietly = TRUE)
  library(ggplot2,    quietly = TRUE)
  library(dplyr,      quietly = TRUE)
  library(tidyr,      quietly = TRUE)
  library(gridExtra,  quietly = TRUE)
  library(grid,       quietly = TRUE)
  library(openxlsx,   quietly = TRUE)
  library(ggnewscale, quietly = TRUE)
  library(data.table, quietly = TRUE)
  library(scales,     quietly = TRUE)
  library(survey,     quietly = TRUE)
  library(robustbase, quietly = TRUE)
  library(cowplot,    quietly = TRUE)
  library(svglite,    quietly = TRUE)
})

# ────────────────────────────────────────────────────────────────
# Nature SI 10.1 + 10.3 一次性 vs μ=50 的统计 helper（v3.5 新增）
# ────────────────────────────────────────────────────────────────

# Sign-flip 置换检验 vs μ=50（Nature SI 10.3 单样本版）
.perm_test_vs_50 <- function(x, n_perm = 10000, mu0 = 50) {
  x <- x[!is.na(x)]
  if (length(x) < 5) return(list(obs = NA_real_, p = NA_real_))
  obs      <- mean(x) - mu0
  centered <- x - mu0
  perm_obs <- replicate(n_perm,
                        mean(centered * sample(c(-1, 1), length(x), replace = TRUE)))
  p <- mean(abs(perm_obs) >= abs(obs))
  if (p == 0) p <- 1 / (n_perm + 1)
  list(obs = obs, p = p)
}

# Cohen's d + Hedges' g (小样本校正) + 95% CI，单样本 vs μ=50
.cohens_d_vs_50 <- function(x, mu0 = 50) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 5) return(list(d = NA, g = NA, ci_lo = NA, ci_hi = NA, se = NA))
  d  <- (mean(x) - mu0) / sd(x)
  J  <- 1 - 3 / (4 * (n - 1) - 1)        # Hedges 小样本校正因子
  g  <- J * d
  se <- sqrt(1/n + g^2 / (2 * n))
  list(d = d, g = g, ci_lo = g - 1.96 * se, ci_hi = g + 1.96 * se, se = se)
}

# Nature SI 10.1 星号映射（基于 BH-FDR P<0.001 等阈值）
.sig_star <- function(p_fdr) {
  dplyr::case_when(
    is.na(p_fdr)  ~ "",
    p_fdr < 0.001 ~ "***",
    p_fdr < 0.01  ~ "**",
    p_fdr < 0.05  ~ "*",
    TRUE          ~ "ns"
  )
}

options(scipen = 999)
options(mc_doScale_quiet = TRUE)   # ← 新增这一行
set.seed(42)   # 保证置换检验可复现

# ========================================================================================
# 0. 全局配置
# ========================================================================================
BASE_DIR     <- "E:/disc/group_meeting4_202602-202608/3/"
OUTPUT_DIR   <- paste0(BASE_DIR, "Part3_Caliber_Results/")
HEALTHY_PATH <- paste0(BASE_DIR, "all_test_with_diameter_top_MERGED.csv")   # ★ 数据源更新
RDS_DIR      <- paste0(OUTPUT_DIR, "RDS_Models/")

# ★★★ 运行模式："model_only" / "disease_analysis" / "plot_only"
RUN_MODE <- "plot_only"

# ════════════════════════════════════════════════════════════════
# 论文级关键图导出：plot_only 时为 4 张关键图额外存矢量 SVG
# ════════════════════════════════════════════════════════════════
SAVE_KEY_SVG <- TRUE   # plot_only 时关键图额外存 SVG

save_key_figure <- function(plot_obj, path_noext, width, height) {
  ggsave(paste0(path_noext, ".pdf"), plot = plot_obj,
         width = width, height = height, dpi = 300,
         device = cairo_pdf, limitsize = FALSE)
  cat(sprintf("  ✅ %s.pdf\n", basename(path_noext)))
  if (isTRUE(SAVE_KEY_SVG) && RUN_MODE == "plot_only") {
    ok <- tryCatch({
      ggsave(paste0(path_noext, ".svg"), plot = plot_obj,
             width = width, height = height,
             device = svglite::svglite, limitsize = FALSE)
      TRUE
    }, error = function(e) {
      cat(sprintf("  ⚠ SVG 保存失败: %s\n", e$message)); FALSE })
    if (ok) cat(sprintf("  ✅ %s.svg  ← 论文用矢量图\n", basename(path_noext)))
  }
}

# ★★★ 多疾病配置
DISEASES <- list(
  list(
    name      = "CKD",
    data_path = "E:/disc/group_meeting4_202602-202608/3/CKD/特征导出/CKD_all_with_diameter_top_MERGED.csv",
    color_female = "#1565C0", color_male = "#64B5F6"
  ),
  list(
    name      = "AD",
    data_path = "E:/disc/group_meeting4_202602-202608/3/CKD/特征导出/AD_all_with_diameter_top_MERGED.csv",
    color_female = "#B71C1C", color_male = "#EF9A9A"
  )
)

SELECTED_SUMMARY_NCOL <- 6

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(RDS_DIR))    dir.create(RDS_DIR,    recursive = TRUE)

cat("\n", strrep("=", 80), "\n")
cat("  Part 3: Caliber（口径，41参数）GAMLSS分析  v15.1\n")
cat(sprintf("  运行模式: %s\n", RUN_MODE))
cat(sprintf("  疾病列表: %s\n", paste(sapply(DISEASES,`[[`,"name"), collapse=", ")))
cat(strrep("=", 80), "\n\n")

# ========================================================================================
# 1. 分析目标（新短标签）
# ========================================================================================

# ───────── 11-12 级（仅 model_only 用；新版主表已不含 L11-12，可整段删除） ─────────
# EXTRA_11_12 <- list(
  # list(col_name="artery11_avg_diameter_um", label="Cal_A_L11", unit="μm", file_tag="artery11_avg_diam", y_group="diameter"),
  # list(col_name="artery12_avg_diameter_um", label="Cal_A_L12", unit="μm", file_tag="artery12_avg_diam", y_group="diameter"),
  # list(col_name="vein11_avg_diameter_um", label="Cal_V_L11", unit="μm", file_tag="vein11_avg_diam", y_group="diameter"),
  # list(col_name="vein12_avg_diameter_um", label="Cal_V_L12", unit="μm", file_tag="vein12_avg_diam", y_group="diameter")
# )

# ───────── Section 1：动脉各分级平均管径（10） ─────────
SEC1_TARGETS <- list(
  list(col_name="artery1_avg_diameter_um", label="Cal_A_L1", unit="μm", file_tag="artery1_avg_diam", y_group="diameter"),
  list(col_name="artery2_avg_diameter_um", label="Cal_A_L2", unit="μm", file_tag="artery2_avg_diam", y_group="diameter"),
  list(col_name="artery3_avg_diameter_um", label="Cal_A_L3", unit="μm", file_tag="artery3_avg_diam", y_group="diameter"),
  list(col_name="artery4_avg_diameter_um", label="Cal_A_L4", unit="μm", file_tag="artery4_avg_diam", y_group="diameter"),
  list(col_name="artery5_avg_diameter_um", label="Cal_A_L5", unit="μm", file_tag="artery5_avg_diam", y_group="diameter"),
  list(col_name="artery6_avg_diameter_um", label="Cal_A_L6", unit="μm", file_tag="artery6_avg_diam", y_group="diameter"),
  list(col_name="artery7_avg_diameter_um", label="Cal_A_L7", unit="μm", file_tag="artery7_avg_diam", y_group="diameter"),
  list(col_name="artery8_avg_diameter_um", label="Cal_A_L8", unit="μm", file_tag="artery8_avg_diam", y_group="diameter"),
  list(col_name="artery9_avg_diameter_um", label="Cal_A_L9", unit="μm", file_tag="artery9_avg_diam", y_group="diameter"),
  list(col_name="artery10_avg_diameter_um", label="Cal_A_L10", unit="μm", file_tag="artery10_avg_diam", y_group="diameter")
)

# ───────── Section 2：静脉各分级平均管径（10） ─────────
SEC2_TARGETS <- list(
  list(col_name="vein1_avg_diameter_um", label="Cal_V_L1", unit="μm", file_tag="vein1_avg_diam", y_group="diameter"),
  list(col_name="vein2_avg_diameter_um", label="Cal_V_L2", unit="μm", file_tag="vein2_avg_diam", y_group="diameter"),
  list(col_name="vein3_avg_diameter_um", label="Cal_V_L3", unit="μm", file_tag="vein3_avg_diam", y_group="diameter"),
  list(col_name="vein4_avg_diameter_um", label="Cal_V_L4", unit="μm", file_tag="vein4_avg_diam", y_group="diameter"),
  list(col_name="vein5_avg_diameter_um", label="Cal_V_L5", unit="μm", file_tag="vein5_avg_diam", y_group="diameter"),
  list(col_name="vein6_avg_diameter_um", label="Cal_V_L6", unit="μm", file_tag="vein6_avg_diam", y_group="diameter"),
  list(col_name="vein7_avg_diameter_um", label="Cal_V_L7", unit="μm", file_tag="vein7_avg_diam", y_group="diameter"),
  list(col_name="vein8_avg_diameter_um", label="Cal_V_L8", unit="μm", file_tag="vein8_avg_diam", y_group="diameter"),
  list(col_name="vein9_avg_diameter_um", label="Cal_V_L9", unit="μm", file_tag="vein9_avg_diam", y_group="diameter"),
  list(col_name="vein10_avg_diameter_um", label="Cal_V_L10", unit="μm", file_tag="vein10_avg_diam", y_group="diameter")
)

# ───────── Section 3：动脉 1-10 级最大血管直径 Top1（10） ─────────
SEC3_TARGETS <- list(
  list(col_name="artery1_diameter_top1_um", label="Cal_A_L1_Top1", unit="μm", file_tag="artery1_diam_top1", y_group="diameter"),
  list(col_name="artery2_diameter_top1_um", label="Cal_A_L2_Top1", unit="μm", file_tag="artery2_diam_top1", y_group="diameter"),
  list(col_name="artery3_diameter_top1_um", label="Cal_A_L3_Top1", unit="μm", file_tag="artery3_diam_top1", y_group="diameter"),
  list(col_name="artery4_diameter_top1_um", label="Cal_A_L4_Top1", unit="μm", file_tag="artery4_diam_top1", y_group="diameter"),
  list(col_name="artery5_diameter_top1_um", label="Cal_A_L5_Top1", unit="μm", file_tag="artery5_diam_top1", y_group="diameter"),
  list(col_name="artery6_diameter_top1_um", label="Cal_A_L6_Top1", unit="μm", file_tag="artery6_diam_top1", y_group="diameter"),
  list(col_name="artery7_diameter_top1_um", label="Cal_A_L7_Top1", unit="μm", file_tag="artery7_diam_top1", y_group="diameter"),
  list(col_name="artery8_diameter_top1_um", label="Cal_A_L8_Top1", unit="μm", file_tag="artery8_diam_top1", y_group="diameter"),
  list(col_name="artery9_diameter_top1_um", label="Cal_A_L9_Top1", unit="μm", file_tag="artery9_diam_top1", y_group="diameter"),
  list(col_name="artery10_diameter_top1_um", label="Cal_A_L10_Top1", unit="μm", file_tag="artery10_diam_top1", y_group="diameter")
)

# ───────── Section 4：静脉 1-10 级最大血管直径 Top1（10） ─────────
SEC4_TARGETS <- list(
  list(col_name="vein1_diameter_top1_um", label="Cal_V_L1_Top1", unit="μm", file_tag="vein1_diam_top1", y_group="diameter"),
  list(col_name="vein2_diameter_top1_um", label="Cal_V_L2_Top1", unit="μm", file_tag="vein2_diam_top1", y_group="diameter"),
  list(col_name="vein3_diameter_top1_um", label="Cal_V_L3_Top1", unit="μm", file_tag="vein3_diam_top1", y_group="diameter"),
  list(col_name="vein4_diameter_top1_um", label="Cal_V_L4_Top1", unit="μm", file_tag="vein4_diam_top1", y_group="diameter"),
  list(col_name="vein5_diameter_top1_um", label="Cal_V_L5_Top1", unit="μm", file_tag="vein5_diam_top1", y_group="diameter"),
  list(col_name="vein6_diameter_top1_um", label="Cal_V_L6_Top1", unit="μm", file_tag="vein6_diam_top1", y_group="diameter"),
  list(col_name="vein7_diameter_top1_um", label="Cal_V_L7_Top1", unit="μm", file_tag="vein7_diam_top1", y_group="diameter"),
  list(col_name="vein8_diameter_top1_um", label="Cal_V_L8_Top1", unit="μm", file_tag="vein8_diam_top1", y_group="diameter"),
  list(col_name="vein9_diameter_top1_um", label="Cal_V_L9_Top1", unit="μm", file_tag="vein9_diam_top1", y_group="diameter"),
  list(col_name="vein10_diameter_top1_um", label="Cal_V_L10_Top1", unit="μm", file_tag="vein10_diam_top1", y_group="diameter")
)

# ───────── Section 5：动脉 1-10 级前6根平均直径 Top6（10） ─────────
SEC5_TARGETS <- list(
  list(col_name="artery1_diameter_top6_avg_um", label="Cal_A_L1_Top6", unit="μm", file_tag="artery1_diam_top6", y_group="diameter"),
  list(col_name="artery2_diameter_top6_avg_um", label="Cal_A_L2_Top6", unit="μm", file_tag="artery2_diam_top6", y_group="diameter"),
  list(col_name="artery3_diameter_top6_avg_um", label="Cal_A_L3_Top6", unit="μm", file_tag="artery3_diam_top6", y_group="diameter"),
  list(col_name="artery4_diameter_top6_avg_um", label="Cal_A_L4_Top6", unit="μm", file_tag="artery4_diam_top6", y_group="diameter"),
  list(col_name="artery5_diameter_top6_avg_um", label="Cal_A_L5_Top6", unit="μm", file_tag="artery5_diam_top6", y_group="diameter"),
  list(col_name="artery6_diameter_top6_avg_um", label="Cal_A_L6_Top6", unit="μm", file_tag="artery6_diam_top6", y_group="diameter"),
  list(col_name="artery7_diameter_top6_avg_um", label="Cal_A_L7_Top6", unit="μm", file_tag="artery7_diam_top6", y_group="diameter"),
  list(col_name="artery8_diameter_top6_avg_um", label="Cal_A_L8_Top6", unit="μm", file_tag="artery8_diam_top6", y_group="diameter"),
  list(col_name="artery9_diameter_top6_avg_um", label="Cal_A_L9_Top6", unit="μm", file_tag="artery9_diam_top6", y_group="diameter"),
  list(col_name="artery10_diameter_top6_avg_um", label="Cal_A_L10_Top6", unit="μm", file_tag="artery10_diam_top6", y_group="diameter")
)

# ───────── Section 6：静脉 1-10 级前6根平均直径 Top6（10） ─────────
SEC6_TARGETS <- list(
  list(col_name="vein1_diameter_top6_avg_um", label="Cal_V_L1_Top6", unit="μm", file_tag="vein1_diam_top6", y_group="diameter"),
  list(col_name="vein2_diameter_top6_avg_um", label="Cal_V_L2_Top6", unit="μm", file_tag="vein2_diam_top6", y_group="diameter"),
  list(col_name="vein3_diameter_top6_avg_um", label="Cal_V_L3_Top6", unit="μm", file_tag="vein3_diam_top6", y_group="diameter"),
  list(col_name="vein4_diameter_top6_avg_um", label="Cal_V_L4_Top6", unit="μm", file_tag="vein4_diam_top6", y_group="diameter"),
  list(col_name="vein5_diameter_top6_avg_um", label="Cal_V_L5_Top6", unit="μm", file_tag="vein5_diam_top6", y_group="diameter"),
  list(col_name="vein6_diameter_top6_avg_um", label="Cal_V_L6_Top6", unit="μm", file_tag="vein6_diam_top6", y_group="diameter"),
  list(col_name="vein7_diameter_top6_avg_um", label="Cal_V_L7_Top6", unit="μm", file_tag="vein7_diam_top6", y_group="diameter"),
  list(col_name="vein8_diameter_top6_avg_um", label="Cal_V_L8_Top6", unit="μm", file_tag="vein8_diam_top6", y_group="diameter"),
  list(col_name="vein9_diameter_top6_avg_um", label="Cal_V_L9_Top6", unit="μm", file_tag="vein9_diam_top6", y_group="diameter"),
  list(col_name="vein10_diameter_top6_avg_um", label="Cal_V_L10_Top6", unit="μm", file_tag="vein10_diam_top6", y_group="diameter")
)

# ───────── Section 7：AVR（6） ─────────
SEC7_TARGETS <- list(
  list(col_name="AVR_L1", label="AVR_L1", unit="", file_tag="AVR_L1", y_group="AVR"),
  list(col_name="AVR_L2", label="AVR_L2", unit="", file_tag="AVR_L2", y_group="AVR"),
  list(col_name="AVR_L3", label="AVR_L3", unit="", file_tag="AVR_L3", y_group="AVR"),
  list(col_name="AVR_L4", label="AVR_L4", unit="", file_tag="AVR_L4", y_group="AVR"),
  list(col_name="AVR_L5", label="AVR_L5", unit="", file_tag="AVR_L5", y_group="AVR"),
  list(col_name="AVR_L6", label="AVR_L6", unit="", file_tag="AVR_L6", y_group="AVR")
)

TARGETS      <- c(SEC1_TARGETS, SEC2_TARGETS, SEC3_TARGETS,
                  SEC4_TARGETS, SEC5_TARGETS, SEC6_TARGETS, SEC7_TARGETS)
TARGETS_FULL <- TARGETS

# ───────── 6 列分段布局（7 节；selected=NULL 表示该节画全部） ─────────
SECTIONS <- list(
  list(label = "Artery Average Diameter per Level (10)",
       tags  = sapply(SEC1_TARGETS, `[[`, "file_tag"), selected = NULL),
  list(label = "Vein Average Diameter per Level (10)",
       tags  = sapply(SEC2_TARGETS, `[[`, "file_tag"), selected = NULL),
  list(label = "Artery Max Diameter Top1 per Level (10)",
       tags  = sapply(SEC3_TARGETS, `[[`, "file_tag"), selected = NULL),
  list(label = "Vein Max Diameter Top1 per Level (10)",
       tags  = sapply(SEC4_TARGETS, `[[`, "file_tag"), selected = NULL),
  list(label = "Artery Top6 Avg Diameter per Level (10)",
       tags  = sapply(SEC5_TARGETS, `[[`, "file_tag"), selected = NULL),
  list(label = "Vein Top6 Avg Diameter per Level (10)",
       tags  = sapply(SEC6_TARGETS, `[[`, "file_tag"), selected = NULL),
  list(label = "AVR (6)",
       tags  = sapply(SEC7_TARGETS, `[[`, "file_tag"), selected = NULL)
)

# ───────── 6 列分段布局（5 节, 30 参数；只画图片展示的子集） ─────────
#SECTIONS <- list(
  #list(label = "Artery Average Diameter (L1–L6)",
       #tags  = sapply(SEC1_TARGETS, `[[`, "file_tag"),
       #selected = c("artery1_avg_diam","artery2_avg_diam","artery3_avg_diam",
                    #"artery4_avg_diam","artery5_avg_diam","artery6_avg_diam")),
  #list(label = "Vein Average Diameter (L1–L6)",
       #tags  = sapply(SEC2_TARGETS, `[[`, "file_tag"),
       #selected = c("vein1_avg_diam","vein2_avg_diam","vein3_avg_diam",
                    #"vein4_avg_diam","vein5_avg_diam","vein6_avg_diam")),
  #list(label = "Artery Max Diameter Top1 (L1–L6)",
       #tags  = sapply(SEC3_TARGETS, `[[`, "file_tag"),
       #selected = c("artery1_diam_top1","artery2_diam_top1","artery3_diam_top1",
                    #"artery4_diam_top1","artery5_diam_top1","artery6_diam_top1")),
  #list(label = "Vein Max Diameter Top1 (L1–L6)",
       #tags  = sapply(SEC4_TARGETS, `[[`, "file_tag"),
       #selected = c("vein1_diam_top1","vein2_diam_top1","vein3_diam_top1",
                    #"vein4_diam_top1","vein5_diam_top1","vein6_diam_top1")),
  #list(label = "AVR (6)",
       #tags  = sapply(SEC7_TARGETS, `[[`, "file_tag"), selected = NULL)
#)

get_section_selected_tags <- function(section) {
  if (is.null(section$selected) || length(section$selected) == 0) section$tags
  else intersect(section$selected, section$tags)
}

cat(sprintf("  TARGETS (用于打分/画图): %d 个，分 %d 节\n",
            length(TARGETS), length(SECTIONS)))
cat(sprintf("  TARGETS_FULL (含11-12级,仅model_only用): %d 个\n\n", length(TARGETS_FULL)))

# ========================================================================================
# 2. 工具函数
# ========================================================================================
FP.SET <- matrix(c(
  1,1,0, 2,1,0, 2,0,1, 2,1,1,
  2,2,1, 2,2,2, 3,1,0, 3,1,1,
  3,2,1, 3,1,2, 3,2,2, 3,3,0,
  3,3,1, 3,3,2
), byrow=TRUE, ncol=3, dimnames=list(NULL, c("mu","sigma","nu")))

# ★★★ Force FO110：仅对下列 file_tag 强制 FO110，其余正常 BIC 选择 ★★★
FORCE_FO110_TAGS <- c("artery1_avg_diam")   # 即 Cal_A_L1
FP.FORCE_110 <- matrix(c(1, 1, 0), byrow = TRUE, ncol = 3,
                       dimnames = list(NULL, c("mu", "sigma", "nu")))

GAMLSS.CTRL <- gamlss.control(n.cyc=400, c.crit=1e-3, trace=FALSE)

find_column <- function(df, possible_names) {
  for (nm in possible_names) if (nm %in% names(df)) return(nm)
  for (nm in possible_names) {
    idx <- which(tolower(names(df)) == tolower(nm))
    if (length(idx) > 0) return(names(df)[idx[1]])
  }
  for (nm in possible_names) {
    idx <- grep(nm, names(df), ignore.case=TRUE)
    if (length(idx) > 0) return(names(df)[idx[1]])
  }
  return(NULL)
}

GG_variance <- function(mu, sigma, nu) {
  AA <- 2*log(mu) - ((2/nu)*(-1)*log((sigma^2*nu^2))) - 2*lgamma(1/(sigma^2*nu^2))
  ww <- lgamma(1/(sigma^2*nu^2)+2/nu) + lgamma(1/(sigma^2*nu^2))
  uu <- 2*lgamma(1/(sigma^2*nu^2)+1/nu)
  BB <- ww + log((1 - exp(uu - ww)))
  ifelse(nu>0|(nu<0 & sigma^2*abs(nu)<0.5),
         ifelse(is.nan(BB), NA, exp(AA+BB)), Inf)
}

remove_outliers_adjbox <- function(x, coef=1.5) {
  stats_res <- tryCatch(
    robustbase::adjboxStats(x, coef=coef),
    error=function(e) { cat("    ⚠ adjboxStats 失败，回退标准围栏\n"); boxplot.stats(x, coef=coef) }
  )
  list(keep=x>=stats_res$fence[1] & x<=stats_res$fence[2],
       lo=stats_res$fence[1], hi=stats_res$fence[2])
}

# ★ ensure_avr：当数据缺少 AVR_1/AVR_2/AVR_3 列时，从 artery/vein top6 直径计算
#   AVR_k = artery_k_diameter_top6_avg / vein_k_diameter_top6_avg
#   优先使用 _um 后缀；若无则尝试 _mm
ensure_avr <- function(df) {
  for (k in 1:6) {                                   # ★ 1:3 → 1:6
    avr_name <- paste0("AVR_L", k)                    # ★ 计算/输出列名统一为 AVR_Lk
    legacy   <- paste0("AVR_", k)                     # 兼容旧命名 AVR_1 / AVR_2 ...
    if (avr_name %in% names(df)) {
      df[[avr_name]] <- suppressWarnings(as.numeric(df[[avr_name]])); next
    }
    if (legacy %in% names(df)) {                      # CSV 若用旧名，则沿用其值
      df[[avr_name]] <- suppressWarnings(as.numeric(df[[legacy]])); next
    }
    a_um <- paste0("artery", k, "_diameter_top6_avg_um")
    v_um <- paste0("vein",   k, "_diameter_top6_avg_um")
    a_mm <- paste0("artery", k, "_diameter_top6_avg_mm")
    v_mm <- paste0("vein",   k, "_diameter_top6_avg_mm")
    if (a_um %in% names(df) && v_um %in% names(df)) {
      a <- suppressWarnings(as.numeric(df[[a_um]]))
      v <- suppressWarnings(as.numeric(df[[v_um]]))
      df[[avr_name]] <- ifelse(is.na(v) | v == 0, NA_real_, a / v)
    } else if (a_mm %in% names(df) && v_mm %in% names(df)) {
      a <- suppressWarnings(as.numeric(df[[a_mm]]))
      v <- suppressWarnings(as.numeric(df[[v_mm]]))
      df[[avr_name]] <- ifelse(is.na(v) | v == 0, NA_real_, a / v)
    } else {
      cat(sprintf("    ⚠ 无法构建 %s（缺少 artery%d/vein%d top6 直径列）\n", avr_name, k, k))
      df[[avr_name]] <- NA_real_
    }
  }
  df
}

fit_gamlss_unified <- function(data, fp_set=FP.SET, ctrl=GAMLSS.CTRL, verbose=TRUE) {
  if (nrow(data) < 20) { cat("    ⚠ 样本量不足，跳过\n"); return(NULL) }
  best_model <- NULL; best_bic <- Inf; best_spec <- NULL; results <- data.frame()
  for (i in 1:nrow(fp_set)) {
    p_mu <- fp_set[i,"mu"]; p_sig <- fp_set[i,"sigma"]; p_nu <- fp_set[i,"nu"]
    spec_str <- sprintf("FO%d%d%d",p_mu,p_sig,p_nu)
    fm <- if(p_mu >0) as.formula(sprintf("Outcome~1+Sex+fp(Age,%d)",p_mu)) else Outcome~1+Sex
    fs <- if(p_sig>0) as.formula(sprintf("~1+Sex+fp(Age,%d)",p_sig))       else ~1+Sex
    fn <- if(p_nu >0) as.formula(sprintf("~1+fp(Age,%d)",p_nu))            else ~1
    m <- tryCatch(
      gamlss(fm,sigma.formula=fs,nu.formula=fn,data=data,
             family=GG(mu.link="log",sigma.link="log",nu.link="identity"),control=ctrl),
      error=function(e){if(verbose)cat(sprintf("    [%s] 失败: %s\n",spec_str,e$message));NULL},
      warning=function(w) suppressWarnings(
        gamlss(fm,sigma.formula=fs,nu.formula=fn,data=data,
               family=GG(mu.link="log",sigma.link="log",nu.link="identity"),control=ctrl))
    )
    if (!is.null(m) && m$converged) {
      bv <- BIC(m); av <- AIC(m)
      if(verbose) cat(sprintf("    [%s] BIC=%9.2f  AIC=%9.2f%s\n",spec_str,bv,av,
                              if(bv<best_bic)"  ←" else ""))
      results <- rbind(results,data.frame(Spec=spec_str,p_mu=as.integer(p_mu),
                                          p_sig=as.integer(p_sig),p_nu=as.integer(p_nu),
                                          BIC=bv,AIC=av,Converged=TRUE))
      if (bv<best_bic){best_bic<-bv;best_model<-m
      best_spec<-list(p_mu=as.integer(p_mu),p_sig=as.integer(p_sig),p_nu=as.integer(p_nu))}
    } else {
      results <- rbind(results,data.frame(Spec=spec_str,p_mu=as.integer(p_mu),
                                          p_sig=as.integer(p_sig),p_nu=as.integer(p_nu),
                                          BIC=NA,AIC=NA,Converged=FALSE))
    }
  }
  if (!is.null(best_model))
    cat(sprintf("  [统一] 最优: FO%d%d%d  BIC=%.2f\n",
                best_spec$p_mu,best_spec$p_sig,best_spec$p_nu,best_bic))
  else cat("  [统一] 所有模型均未收敛\n")
  list(model=best_model,bic=best_bic,spec=best_spec,selection_table=results)
}

predict_centiles_unified <- function(fit_result, age_seq, sex_val, fit_data,
                                     probs=c(0.025,0.25,0.5,0.75,0.975)) {
  model <- fit_result$model
  empty <- data.frame(Age=age_seq,P025=NA,P25=NA,P50=NA,P75=NA,P975=NA,Sex=sex_val)
  if (is.null(model)) return(empty)
  pred <- tryCatch(
    predictAll(model,newdata=data.frame(Age=age_seq,Sex=factor(sex_val,levels=levels(fit_data$Sex))),
               data=fit_data,type="response"),
    error=function(e) NULL)
  if (is.null(pred)) return(empty)
  q_mat <- mapply(function(mu,sigma,nu)
    tryCatch(qGG(probs,mu=mu,sigma=sigma,nu=nu),error=function(e)rep(NA_real_,5)),
    pred$mu,pred$sigma,pred$nu)
  data.frame(Age=age_seq,P025=q_mat[1,],P25=q_mat[2,],P50=q_mat[3,],
             P75=q_mat[4,],P975=q_mat[5,],Sex=sex_val)
}

calculate_variance_unified <- function(fit_result, age_seq, sex_val, fit_data) {
  model <- fit_result$model
  if (is.null(model)) return(data.frame(Age=age_seq,Variance=NA,Sex=sex_val))
  pred <- tryCatch(
    predictAll(model,newdata=data.frame(Age=age_seq,Sex=factor(sex_val,levels=levels(fit_data$Sex))),
               data=fit_data,type="response"),
    error=function(e) NULL)
  if (is.null(pred)) return(data.frame(Age=age_seq,Variance=NA,Sex=sex_val))
  data.frame(Age=age_seq,Variance=mapply(GG_variance,pred$mu,pred$sigma,pred$nu),Sex=sex_val)
}

.calc_velocity_core <- function(predict_fun, age_seq) {
  age_dense <- seq(min(age_seq),max(age_seq),by=0.05)
  med_dense  <- predict_fun(age_dense)
  if (all(is.na(med_dense))) return(rep(NA_real_,length(age_seq)))
  dx <- diff(age_dense); dy <- diff(med_dense); vel <- dy/dx
  vel_full <- c(vel[1],vel)
  approx(x=age_dense,y=vel_full,xout=age_seq,rule=2)$y
}

calculate_velocity_unified <- function(fit_result, age_seq, sex_val, fit_data) {
  model <- fit_result$model
  if (is.null(model)) return(data.frame(Age=age_seq,Velocity=NA,Sex=sex_val))
  predict_fun <- function(age_vec) {
    pred <- tryCatch(
      predictAll(model,newdata=data.frame(Age=age_vec,Sex=factor(sex_val,levels=levels(fit_data$Sex))),
                 data=fit_data,type="response"),
      error=function(e) NULL)
    if (is.null(pred)) return(rep(NA_real_,length(age_vec)))
    mapply(function(mu,sigma,nu)
      tryCatch(qGG(0.5,mu=mu,sigma=sigma,nu=nu),error=function(e) NA_real_),
      pred$mu,pred$sigma,pred$nu)
  }
  data.frame(Age=age_seq,Velocity=.calc_velocity_core(predict_fun,age_seq),Sex=sex_val)
}

extract_polynomial_formula_nature <- function(adaptive_result, var_label="Variable") {
  parse_one <- function(model, spec, sex_note) {
    if (is.null(model)||is.null(spec)) return(NULL)
    result <- list(spec=sprintf("FO%d%d%d",spec$p_mu,spec$p_sig,spec$p_nu),sex=sex_note)
    for (LAB in c("mu","sigma","nu")) {
      b <- tryCatch(coef(model,what=LAB),error=function(e) NULL)
      if (is.null(b)){result[[LAB]] <- paste0(LAB,": (failed)");next}
      lhs <- switch(LAB,mu="log(\u03bc)",sigma="log(\u03c3)",nu="\u03bd")
      fp_powers <- NULL; fp_coefs <- NULL
      tryCatch({
        smo_list <- getSmo(model,what=LAB,which=0)
        if (!is.null(smo_list)&&length(smo_list)>0) {
          CLASS <- sapply(smo_list,class)
          if (any(CLASS=="lm")) {
            fp_smo <- getSmo(model,what=LAB,which=which(CLASS=="lm")[1])
            if (!is.null(fp_smo$power))        fp_powers <- fp_smo$power
            if (!is.null(fp_smo$coefficients)) fp_coefs  <- fp_smo$coefficients
            else {
              age_idx <- setdiff(grep("fp\\(|Age|age",names(b),ignore.case=TRUE),
                                 grep("Intercept|Sex|sex",names(b),ignore.case=TRUE))
              if (length(age_idx)>0) fp_coefs <- b[age_idx]
            }
          }
        }
      }, error=function(e) {
        age_idx <- setdiff(grep("fp\\(|Age|age",names(b),ignore.case=TRUE),
                           grep("Intercept|Sex|sex",names(b),ignore.case=TRUE))
        if (length(age_idx)>0) fp_coefs <<- b[age_idx]
      })
      terms <- character(0)
      if ("(Intercept)" %in% names(b)) terms <- c(terms,sprintf("\u03b1_%s",LAB))
      sex_idx <- grep("^Sex",names(b))
      if (length(sex_idx)>0) terms <- c(terms,sprintf("\u03b1_{%s,sex}(sex)",LAB))
      if (!is.null(fp_powers)&&!is.null(fp_coefs)&&length(fp_coefs)>0) {
        n_t <- min(length(fp_powers),length(fp_coefs))
        for (k in seq_len(n_t)) {
          pw <- fp_powers[k]; is_rep <- if(k>1) pw==fp_powers[k-1] else FALSE
          pw_str <- if(pw==round(pw)){
            pi <- as.integer(pw)
            if(pi==1)"(age)" else if(pi==0)"log(age)" else
              if(pi==-1)"(age)^{-1}" else sprintf("(age)^{%d}",pi)
          } else sprintf("(age)^{%.2f}",pw)
          if (is_rep) pw_str <- paste0(pw_str,"\u00b7log(age)")
          terms <- c(terms,sprintf("\u03b2_{%s,%d}%s",LAB,k,pw_str))
        }
      } else if (length(fp_coefs)>0&&is.null(fp_powers)) {
        for (k in seq_along(fp_coefs))
          terms <- c(terms,sprintf("\u03b2_{%s,%d}(age)^{?}",LAB,k))
      }
      result[[LAB]] <- sprintf("%s = %s",lhs,paste(terms,collapse=" + "))
      coef_details <- character(0)
      if ("(Intercept)" %in% names(b))
        coef_details <- c(coef_details,sprintf("  \u03b1_%s = %.6f",LAB,b["(Intercept)"]))
      if (length(sex_idx)>0)
        coef_details <- c(coef_details,sprintf("  \u03b1_{%s,sex} = %.6f",LAB,b[sex_idx[1]]))
      if (!is.null(fp_coefs)&&length(fp_coefs)>0)
        for (k in seq_along(fp_coefs))
          coef_details <- c(coef_details,
                            sprintf("  \u03b2_{%s,%d} = %s",LAB,k,
                                    if(is.na(fp_coefs[k]))"NA" else sprintf("%.8f",fp_coefs[k])))
      result[[paste0(LAB,"_coefficients")]] <- paste(coef_details,collapse="\n")
    }
    return(result)
  }
  list(type="unified",var_label=var_label,
       unified=parse_one(adaptive_result$unified$model,adaptive_result$unified$spec,"Unified model (both sexes)"))
}

write_polynomial_formula_nature <- function(formula_output, label, file_tag, output_dir) {
  out_path <- paste0(output_dir,"Formula_Nature_",file_tag,".txt")
  sink(file=out_path)
  cat(strrep("=",80),"\n")
  cat(sprintf("Variable: %s\n",label)); cat(sprintf("Model Type: %s\n",formula_output$type))
  cat("Format: Nature publication style\n"); cat(strrep("=",80),"\n\n")
  cat(sprintf("%s ~ Generalised Gamma(\u03bc, \u03c3, \u03bd) with\n\n",label))
  print_one <- function(res,title) {
    if (is.null(res)){cat(sprintf("\n%s: Model did not converge\n\n",title));return(invisible())}
    cat(sprintf("%s\n",title)); cat(strrep("-",70),"\n")
    cat(sprintf("Specification: %s\n",res$spec)); cat(sprintf("Sex group: %s\n\n",res$sex))
    for (LAB in c("mu","sigma","nu")) if(!is.null(res[[LAB]])) cat(sprintf("%s\n\n",res[[LAB]]))
    cat("Coefficient values:\n"); cat(strrep("-",40),"\n")
    for (LAB in c("mu","sigma","nu")) {
      ck <- paste0(LAB,"_coefficients")
      if(!is.null(res[[ck]])&&nzchar(res[[ck]])) cat(sprintf("%s:\n%s\n\n",toupper(LAB),res[[ck]]))
    }
  }
  print_one(formula_output$unified,"UNIFIED MODEL")
  sink()
  cat(sprintf("  ✅ Formula_Nature_%s.txt\n",file_tag))
}

.score_from_pred <- function(value, pred) {
  mu <- pred$mu; sigma <- pred$sigma; nu <- pred$nu
  pct  <- tryCatch(pGG(value,mu=mu,sigma=sigma,nu=nu)*100, error=function(e) NA)
  med  <- tryCatch(qGG(0.5,mu=mu,sigma=sigma,nu=nu),       error=function(e) NA)
  sd_v <- sqrt(GG_variance(mu,sigma,nu)); zsc <- (value-med)/sd_v
  cat_v <- dplyr::case_when(
    is.na(pct) ~ "无法计算",   pct<2.5  ~ "显著偏低 (<P2.5)",
    pct>97.5   ~ "显著偏高 (>P97.5)", pct<25 ~ "偏低 (<P25)",
    pct>75     ~ "偏高 (>P75)",       TRUE   ~ "正常 (P25-P75)")
  list(percentile=round(pct,2),zscore=round(zsc,3),category=cat_v,
       predicted_median=round(med,4),predicted_sd=round(sd_v,4))
}

calculate_percentile_adaptive <- function(value, age, sex, adaptive_result, df_h) {
  model <- adaptive_result$unified$model
  if (is.null(model)) return(list(percentile=NA,zscore=NA,category="模型不可用",
                                  predicted_median=NA,predicted_sd=NA))
  pred <- tryCatch(
    predictAll(model,newdata=data.frame(Age=age,Sex=factor(sex,levels=levels(df_h$Sex))),
               data=df_h,type="response"),
    error=function(e) NULL)
  if (is.null(pred)) return(list(percentile=NA,zscore=NA,category="预测失败",
                                 predicted_median=NA,predicted_sd=NA))
  .score_from_pred(value,pred)
}

# ========================================================================================
# 3. 阶段1：建模阶段
# ========================================================================================
ALL_RESULTS <- list()
outlier_log <- list()
model_selection_log <- list()

if (RUN_MODE == "model_only") {
  cat("\n>>> 阶段1: GAMLSS 建模 + 健康参考输出 <<<\n")
  
  cat("--- 读取健康人群数据 ---\n")
  df_healthy_raw <- if (grepl("\\.csv$", HEALTHY_PATH, ignore.case=TRUE)) {
    as.data.frame(data.table::fread(HEALTHY_PATH, encoding="UTF-8"))
  } else {
    openxlsx::read.xlsx(HEALTHY_PATH)
  }
  cat(sprintf("健康人群: %d 行 × %d 列\n", nrow(df_healthy_raw), ncol(df_healthy_raw)))
  
  # ★ 健康人群：补 AVR_1/2/3
  df_healthy_raw <- ensure_avr(df_healthy_raw)
  
  for (i in seq_along(TARGETS_FULL)) {
    tgt      <- TARGETS_FULL[[i]]
    col_name <- tgt$col_name; label    <- tgt$label
    file_tag <- tgt$file_tag; unit     <- tgt$unit
    
    cat("\n",strrep("=",80),"\n")
    cat(sprintf("  [%d/%d] %s\n",i,length(TARGETS_FULL),label))
    cat(strrep("=",80),"\n")
    
    if (!col_name %in% colnames(df_healthy_raw)) {
      cat(sprintf("  ⚠ 健康数据缺少列 '%s'，跳过\n",col_name)); next
    }
    
    df_h_raw <- df_healthy_raw %>%
      select(all_of(c("Image_name",col_name,"Sex","Age"))) %>%
      rename(Outcome=all_of(col_name)) %>%
      mutate(Outcome=as.numeric(as.character(Outcome)),
             Age=as.numeric(as.character(Age)),
             Sex=dplyr::case_when(
               tolower(as.character(Sex)) %in% c("female","女","f","2") ~ "Female",
               tolower(as.character(Sex)) %in% c("male","男","m","1")   ~ "Male",
               TRUE ~ NA_character_)) %>%
      drop_na(Outcome,Sex,Age) %>%
      mutate(Sex=factor(Sex,levels=c("Female","Male")))
    
    n_raw_total <- nrow(df_h_raw); n_raw_f <- sum(df_h_raw$Sex=="Female"); n_raw_m <- sum(df_h_raw$Sex=="Male")
    keep_f <- if(n_raw_f>=10) remove_outliers_adjbox(df_h_raw$Outcome[df_h_raw$Sex=="Female"])$keep else rep(TRUE,n_raw_f)
    keep_m <- if(n_raw_m>=10) remove_outliers_adjbox(df_h_raw$Outcome[df_h_raw$Sex=="Male"])$keep   else rep(TRUE,n_raw_m)
    keep_vec <- logical(n_raw_total)
    keep_vec[df_h_raw$Sex=="Female"] <- keep_f; keep_vec[df_h_raw$Sex=="Male"] <- keep_m
    df_h <- df_h_raw[keep_vec,]
    n_clean_total <- nrow(df_h); n_clean_f <- sum(df_h$Sex=="Female"); n_clean_m <- sum(df_h$Sex=="Male")
    
    outlier_log[[file_tag]] <- data.frame(
      Variable=label,ColName=col_name,
      Outlier_Method="robustbase::adjboxStats(coef=1.5)",
      Raw_Total=n_raw_total,Raw_Female=n_raw_f,Raw_Male=n_raw_m,
      Clean_Total=n_clean_total,Clean_Female=n_clean_f,Clean_Male=n_clean_m,
      Removed_Total=n_raw_total-n_clean_total,
      Removed_Female=n_raw_f-n_clean_f,Removed_Male=n_raw_m-n_clean_m,
      Removed_Pct_Total=round(100*(n_raw_total-n_clean_total)/max(n_raw_total,1),2),
      Removed_Pct_Female=round(100*(n_raw_f-n_clean_f)/max(n_raw_f,1),2),
      Removed_Pct_Male=round(100*(n_raw_m-n_clean_m)/max(n_raw_m,1),2),
      stringsAsFactors=FALSE)
    cat(sprintf("\n  【异常值】原始%d → 保留%d (剔除%d, %.1f%%)\n",
                n_raw_total,n_clean_total,n_raw_total-n_clean_total,
                100*(n_raw_total-n_clean_total)/max(n_raw_total,1)))
    
    # =========================================================================================
    cat("\n--- 统一模型建模 ---\n")
    fp_set_use <- if (file_tag %in% FORCE_FO110_TAGS) {
      cat("  ★ 该参数强制 FO110（mu=fp(Age,1), sigma=fp(Age,1), nu=const）\n")
      FP.FORCE_110
    } else FP.SET
    unified_result <- fit_gamlss_unified(df_h, fp_set_use, GAMLSS.CTRL, verbose=TRUE)
    if (is.null(unified_result$model)) { cat("  ✗ 模型未收敛，跳过\n"); next }
    
    adaptive_result <- list(type="unified",unified=unified_result,female=NULL,male=NULL)
    spec <- unified_result$spec
    model_selection_log[[file_tag]] <- data.frame(
      Variable=label,Model_Type="Unified",BIC=round(unified_result$bic,2),
      FO_Spec=sprintf("FO%d%d%d",spec$p_mu,spec$p_sig,spec$p_nu),stringsAsFactors=FALSE)
    
    write.xlsx(unified_result$selection_table,
               paste0(OUTPUT_DIR,"ModelSel_",file_tag,".xlsx"),rowNames=FALSE)
    cat(sprintf("  ✅ ModelSel_%s.xlsx\n",file_tag))
    
    cat("\n--- 提取多项式公式 ---\n")
    formula_output <- extract_polynomial_formula_nature(adaptive_result,var_label=label)
    write_polynomial_formula_nature(formula_output,label,file_tag,OUTPUT_DIR)
    
    age_min <- quantile(df_h$Age,0.05); age_max <- quantile(df_h$Age,0.95)
    age_seq <- seq(age_min,age_max,by=0.5)
    pred_healthy <- bind_rows(
      predict_centiles_unified(unified_result,age_seq,"Female",df_h),
      predict_centiles_unified(unified_result,age_seq,"Male",  df_h)) %>% drop_na(P50)
    variance_data <- bind_rows(
      calculate_variance_unified(unified_result,age_seq,"Female",df_h),
      calculate_variance_unified(unified_result,age_seq,"Male",  df_h)) %>% drop_na(Variance)
    velocity_data <- bind_rows(
      calculate_velocity_unified(unified_result,age_seq,"Female",df_h),
      calculate_velocity_unified(unified_result,age_seq,"Male",  df_h)) %>% drop_na(Velocity)
    model_note <- sprintf("Unified model  |  FO%d%d%d  |  BIC=%.1f",
                          spec$p_mu,spec$p_sig,spec$p_nu,unified_result$bic)
    plot_pts <- df_h %>% filter(Age>=age_min,Age<=age_max)
    
    cat("\n--- 参考百分位数 Excel ---\n")
    age_ref <- seq(floor(age_min),ceiling(age_max),by=1)
    ref_list <- list()
    for (sx in c("Female","Male")) {
      cr <- predict_centiles_unified(unified_result,age_ref,sx,df_h)
      if (any(!is.na(cr$P50)))
        ref_list[[sx]] <- data.frame(Sex=sx,Age=cr$Age,
                                     P2.5=round(cr$P025,4),P25=round(cr$P25,4),P50=round(cr$P50,4),
                                     P75=round(cr$P75,4),P97.5=round(cr$P975,4),
                                     IQR=round(cr$P75-cr$P25,4),Width95=round(cr$P975-cr$P025,4))
    }
    if (length(ref_list)>0) {
      wb2 <- createWorkbook()
      hs  <- createStyle(fontColour="#FFFFFF",bgFill="#4472C4",halign="center",
                         textDecoration="bold",border="Bottom")
      for (nm in names(ref_list)) {
        addWorksheet(wb2,nm); writeData(wb2,nm,ref_list[[nm]])
        addStyle(wb2,nm,hs,rows=1,cols=1:ncol(ref_list[[nm]]),gridExpand=TRUE)
        setColWidths(wb2,nm,cols=1:ncol(ref_list[[nm]]),widths="auto")
      }
      saveWorkbook(wb2,paste0(OUTPUT_DIR,"Reference_Centiles_",file_tag,".xlsx"),overwrite=TRUE)
      cat(sprintf("  ✅ Reference_Centiles_%s.xlsx\n",file_tag))
    }
    
    result_obj <- list(
      pred_data      = pred_healthy,  variance_data  = variance_data,
      velocity_data  = velocity_data, plot_pts       = plot_pts,
      label          = label,         unit           = unit,
      model_note     = model_note,    adaptive_result= adaptive_result,
      df_h           = df_h)
    ALL_RESULTS[[file_tag]] <- result_obj
    saveRDS(result_obj,paste0(RDS_DIR,"result_",file_tag,".rds"))
    cat(sprintf("  ✅ RDS 已保存: result_%s.rds\n",file_tag))
  }
  
  if (length(outlier_log)>0) {
    outlier_df <- bind_rows(lapply(outlier_log,as.data.frame))
    wb_out <- createWorkbook()
    hs_out <- createStyle(fontColour="#FFFFFF",bgFill="#375623",halign="center",
                          textDecoration="bold",border="Bottom")
    addWorksheet(wb_out,"Outlier_Removal_Summary")
    writeData(wb_out,"Outlier_Removal_Summary",outlier_df)
    addStyle(wb_out,"Outlier_Removal_Summary",hs_out,rows=1,cols=1:ncol(outlier_df),gridExpand=TRUE)
    setColWidths(wb_out,"Outlier_Removal_Summary",cols=1:ncol(outlier_df),widths="auto")
    saveWorkbook(wb_out,paste0(OUTPUT_DIR,"Outlier_Removal_Report_Part3.xlsx"),overwrite=TRUE)
    cat("\n✅ Outlier_Removal_Report_Part3.xlsx\n")
  }
  
  if (length(model_selection_log)>0) {
    model_log_df <- bind_rows(lapply(model_selection_log,as.data.frame))
    write.xlsx(model_log_df, paste0(OUTPUT_DIR,"Model_Selection_Log_Part3.xlsx"), rowNames=FALSE)
    cat("✅ Model_Selection_Log_Part3.xlsx\n")
  }
  
  cat("\n✅ 建模阶段完成！RDS 保存在:", RDS_DIR, "\n")
  
} else {
  cat("\n>>> 从 RDS 加载已保存的健康模型（仅 41 个分析变量）<<<\n")
  for (tgt in TARGETS) {
    rds_path <- paste0(RDS_DIR,"result_",tgt$file_tag,".rds")
    if (file.exists(rds_path)) {
      ALL_RESULTS[[tgt$file_tag]] <- readRDS(rds_path)
    } else {
      cat(sprintf("  ⚠ 未找到RDS: %s（请先运行 model_only 模式）\n",tgt$file_tag))
    }
  }
  cat(sprintf("  ✅ 已加载 %d 个模型\n", length(ALL_RESULTS)))
}

# ★★★ v15.1：用 TARGETS 中定义的新短标签覆盖 RDS 中的旧 label
cat("\n>>> 用新短标签覆盖已加载模型中的旧标签 <<<\n")
override_count <- 0
for (tgt in TARGETS) {
  if (!is.null(ALL_RESULTS[[tgt$file_tag]])) {
    ALL_RESULTS[[tgt$file_tag]]$label <- tgt$label
    override_count <- override_count + 1
  }
}
cat(sprintf("  ✅ 已覆盖 %d 个模型的 label 字段\n", override_count))

# ========================================================================================
# 4. 阶段2：多疾病分析
# ========================================================================================
disease_data_all <- list()

if (RUN_MODE == "disease_analysis") {
  cat("\n>>> 阶段2: 多疾病数据打分 <<<\n")
  
  for (DIS in DISEASES) {
    cat("\n", strrep("─",80), "\n")
    cat(sprintf(">>> 处理疾病: %s\n", DIS$name))
    cat(sprintf("    数据路径: %s\n", DIS$data_path))
    cat(strrep("─",80), "\n")
    
    DISEASE_DIR <- paste0(OUTPUT_DIR, DIS$name, "/")
    if (!dir.exists(DISEASE_DIR)) dir.create(DISEASE_DIR, recursive = TRUE)
    
    if (!file.exists(DIS$data_path)) {
      cat(sprintf("  ⚠ 数据文件不存在，跳过 %s\n", DIS$name)); next
    }
    
    if (grepl("\\.csv$", DIS$data_path, ignore.case=TRUE)) {
      df_dis_raw <- data.table::fread(DIS$data_path, encoding="UTF-8")
      df_dis_raw <- as.data.frame(df_dis_raw)
    } else {
      df_dis_raw <- openxlsx::read.xlsx(DIS$data_path)
    }
    cat(sprintf("  患者数据: %d 行 × %d 列\n", nrow(df_dis_raw), ncol(df_dis_raw)))
    
    # ★ 患者数据：补 AVR_1/2/3
    df_dis_raw <- ensure_avr(df_dis_raw)
    
    age_col <- find_column(df_dis_raw, c("Age","age","年龄"))
    sex_col <- find_column(df_dis_raw, c("Sex","sex","性别","Gender"))
    id_col  <- find_column(df_dis_raw, c("Image_name","PatientID","Patient_ID","ID","id"))
    if (is.null(age_col) || is.null(sex_col)) {
      cat("  ⚠ 缺少 Age 或 Sex 列，跳过\n"); next
    }
    
    all_patient_reports <- list()
    
    for (i in seq_along(TARGETS)) {
      tgt <- TARGETS[[i]]
      col_name <- tgt$col_name; label <- tgt$label; file_tag <- tgt$file_tag
      res <- ALL_RESULTS[[file_tag]]
      if (is.null(res)) { cat(sprintf("  ⚠ [%s] 缺少健康模型\n",label)); next }
      
      cat(sprintf("\n  [%d/%d] %s 患者打分中...\n",i,length(TARGETS),label))
      df_dis <- df_dis_raw %>%
        mutate(
          Outcome   = if(col_name %in% names(df_dis_raw)) as.numeric(get(col_name)) else NA_real_,
          Age       = as.numeric(get(age_col)),
          Sex       = dplyr::case_when(
            tolower(as.character(get(sex_col))) %in% c("female","女","f","2") ~ "Female",
            tolower(as.character(get(sex_col))) %in% c("male","男","m","1")   ~ "Male",
            TRUE ~ NA_character_),
          PatientID = if(!is.null(id_col)) as.character(get(id_col)) else as.character(dplyr::row_number()),
          Disease   = DIS$name
        ) %>% filter(!is.na(Outcome),!is.na(Age),!is.na(Sex)) %>%
        select(PatientID,Age,Sex,Outcome,Disease)
      
      if (nrow(df_dis)==0) { cat("    ⚠ 无有效数据\n"); next }
      
      patient_scores <- df_dis %>% rowwise() %>%
        mutate(score_result=list(calculate_percentile_adaptive(
          Outcome,Age,Sex,res$adaptive_result,res$df_h))) %>%
        ungroup() %>%
        mutate(Percentile      =sapply(score_result,`[[`,"percentile"),
               Zscore          =sapply(score_result,`[[`,"zscore"),
               Category        =sapply(score_result,`[[`,"category"),
               Predicted_Median=sapply(score_result,`[[`,"predicted_median"),
               Predicted_SD    =sapply(score_result,`[[`,"predicted_sd")) %>%
        select(-score_result)
      
      all_patient_reports[[file_tag]] <- patient_scores
      write.xlsx(patient_scores, paste0(DISEASE_DIR,"Patient_Percentile_",file_tag,".xlsx"), rowNames=FALSE)
      cat(sprintf("    ✅ Patient_Percentile_%s.xlsx (n=%d)\n",file_tag,nrow(patient_scores)))
    }
    
    if (length(all_patient_reports)>0) {
      df_all_scores <- bind_rows(lapply(names(all_patient_reports), function(tag)
        all_patient_reports[[tag]] %>%
          select(PatientID,Age,Sex,Disease,Outcome,Percentile,Zscore,Category) %>%
          mutate(Variable=tag)))
      
      patient_summary <- df_all_scores %>%
        group_by(PatientID,Age,Sex,Disease) %>%
        summarise(Total=n(),
                  N_Sig =sum(Category %in% c("显著偏低 (<P2.5)","显著偏高 (>P97.5)")),
                  N_Mild=sum(Category %in% c("偏低 (<P25)","偏高 (>P75)")),
                  N_Norm=sum(Category=="正常 (P25-P75)"),
                  Abn_Rate=round(100*(N_Sig+N_Mild)/Total,1),.groups="drop") %>%
        arrange(desc(N_Sig))
      sex_summary <- df_all_scores %>%
        filter(!is.na(Percentile)) %>%
        group_by(Disease,Variable,Sex) %>%
        summarise(N=n(),Med_Pct=round(median(Percentile,na.rm=TRUE),1),
                  Pct_Low=round(mean(Percentile<2.5,na.rm=TRUE)*100,1),
                  Pct_High=round(mean(Percentile>97.5,na.rm=TRUE)*100,1),.groups="drop")
      
      wb_f <- createWorkbook()
      hs2  <- createStyle(fontColour="#FFFFFF",bgFill="#2E4057",halign="center",
                          textDecoration="bold",border="Bottom")
      add_sh <- function(wb,nm,dat){
        addWorksheet(wb,nm); writeData(wb,nm,dat)
        addStyle(wb,nm,hs2,rows=1,cols=1:ncol(dat),gridExpand=TRUE)
        setColWidths(wb,nm,cols=1:ncol(dat),widths="auto")
      }
      add_sh(wb_f,"Patient_Summary",patient_summary)
      add_sh(wb_f,"All_Scores",     df_all_scores)
      add_sh(wb_f,"Sex_Summary",    sex_summary)
      saveWorkbook(wb_f, paste0(DISEASE_DIR,"Part3_",DIS$name,"_Summary.xlsx"), overwrite=TRUE)
      cat(sprintf("\n  ✅ Part3_%s_Summary.xlsx\n", DIS$name))
      
      saveRDS(list(
        df_all_scores       = df_all_scores,
        all_patient_reports = all_patient_reports,
        disease_name        = DIS$name,
        config              = DIS
      ), paste0(RDS_DIR,"disease_scores_",DIS$name,".rds"))
      cat(sprintf("  ✅ disease_scores_%s.rds\n", DIS$name))
      
      disease_data_all[[DIS$name]] <- list(
        df_all_scores       = df_all_scores,
        all_patient_reports = all_patient_reports,
        config              = DIS
      )
    }
  }
  
} else if (RUN_MODE == "plot_only") {
  cat("\n>>> 从 RDS 加载所有疾病打分数据 <<<\n")
  for (DIS in DISEASES) {
    rds_path <- paste0(RDS_DIR,"disease_scores_",DIS$name,".rds")
    if (file.exists(rds_path)) {
      dis_obj <- readRDS(rds_path)
      disease_data_all[[DIS$name]] <- list(
        df_all_scores       = dis_obj$df_all_scores,
        all_patient_reports = dis_obj$all_patient_reports,
        config              = DIS
      )
      cat(sprintf("  ✅ %s: %d 个变量\n", DIS$name, length(dis_obj$all_patient_reports)))
    } else {
      cat(sprintf("  ⚠ 未找到 %s 的 RDS\n", DIS$name))
    }
  }
}

# ========================================================================================
# 5. 画图函数
# ========================================================================================

# ── 小提琴图（健康人群）
plot_violins <- function(df_model, label, file_tag, unit) {
  y_lab <- if(nchar(unit)>0) paste0(label," (",unit,")") else label
  age_t <- quantile(df_model$Age,c(1/3,2/3))
  df_vio <- df_model %>%
    mutate(AgeGroup=case_when(
      Age<age_t[1]  ~ paste0("<",round(age_t[1]),"yr"),
      Age<=age_t[2] ~ paste0(round(age_t[1]),"-",round(age_t[2]),"yr"),
      TRUE          ~ paste0(">",round(age_t[2]),"yr")),
      AgeGroup=factor(AgeGroup,levels=c(paste0("<",round(age_t[1]),"yr"),
                                        paste0(round(age_t[1]),"-",round(age_t[2]),"yr"),
                                        paste0(">",round(age_t[2]),"yr"))))
  make_annot <- function(dat)
    dat %>% group_by(Sex) %>%
    summarise(Q1=round(quantile(Outcome,.25),3),Med=round(median(Outcome),3),
              Q3=round(quantile(Outcome,.75),3),y_pos=median(Outcome),.groups="drop") %>%
    mutate(label=paste0("Q1:",Q1,"\nMed:",Med,"\nQ3:",Q3))
  get_pval <- function(dat) {
    if(nrow(dat)<10) return("ns")
    r <- tryCatch(svyranktest(Outcome~Sex,design=svydesign(ids=~1,weights=~1,data=dat),
                              test="wilcoxon"),error=function(e)NULL)
    if(is.null(r)) return("ns")
    p <- r$p.value
    if(p<0.001)"***" else if(p<0.01)"**" else if(p<0.05)"*" else "ns"
  }
  one_vio <- function(dat,title_str,show_y=TRUE) {
    ann <- make_annot(dat); ys <- max(dat$Outcome)*1.01
    ggplot(dat,aes(x=Sex,y=Outcome,fill=Sex)) +
      geom_violin(alpha=0.75,trim=TRUE) + geom_boxplot(width=0.12,color="black",outlier.shape=NA) +
      geom_text(data=ann,aes(x=Sex,y=y_pos,label=label),size=2.5,color="white",fontface="bold",lineheight=1.2) +
      annotate("text",x=1.5,y=ys,label=get_pval(dat),size=5,fontface="bold") +
      scale_fill_manual(values=c(Female="#E8595A",Male="#4BACC6")) +
      labs(title=title_str,x=NULL,y=if(show_y)y_lab else NULL) +
      theme_bw() + theme(plot.title=element_text(hjust=0.5,size=11,face="bold"),
                         axis.title.y=element_text(size=9),legend.position="none",panel.grid=element_blank())
  }
  lvls <- levels(df_vio$AgeGroup)
  p_comb <- gridExtra::grid.arrange(
    one_vio(df_vio,"Overall sample",TRUE),
    one_vio(df_vio%>%filter(AgeGroup==lvls[1]),lvls[1],FALSE),
    one_vio(df_vio%>%filter(AgeGroup==lvls[2]),lvls[2],FALSE),
    one_vio(df_vio%>%filter(AgeGroup==lvls[3]),lvls[3],FALSE),
    nrow=1,
    top=grid::textGrob(paste0(label," by Sex and Age Group"),gp=grid::gpar(fontface="bold",fontsize=13)),
    bottom=grid::textGrob(paste0("n=",nrow(df_model),"  |  *** p<0.001  ** p<0.01  * p<0.05"),
                          gp=grid::gpar(fontsize=9,col="#555")))
  ggsave(paste0(OUTPUT_DIR,"ViolinPlot_",file_tag,".pdf"),plot=p_comb,width=16,height=7,dpi=300,device=cairo_pdf)
  cat(sprintf("  ✅ ViolinPlot_%s.pdf\n",file_tag)); return(invisible(p_comb))
}

# ── Trajectories
plot_trajectories <- function(pred_data, plot_pts, label, file_tag, unit, y_lim=NULL, model_note="") {
  if (is.null(pred_data)||nrow(pred_data)==0){cat(sprintf("  ⚠ 跳过 Trajectories [%s]\n",file_tag));return(invisible(NULL))}
  col_f <- "#d62728"; col_m <- "#1f77b4"
  y_lab <- if(nchar(unit)>0) paste0(label," (",unit,")") else label
  # ★ 用 base R 安全取 Female 末端点（替代 slice_max，规避 dplyr 1.1+ 边界 bug）
  .fem <- pred_data[!is.na(pred_data$Age) & pred_data$Sex == "Female", , drop = FALSE]
  female_end <- if (nrow(.fem) > 0) .fem[which.max(.fem$Age), , drop = FALSE] else pred_data[0, , drop = FALSE]
  n_f <- sum(plot_pts$Sex=="Female"); n_m <- sum(plot_pts$Sex=="Male")
  if (is.null(y_lim)) {
    bv <- c(pred_data$P25,pred_data$P75); bv <- bv[!is.na(bv)]
    if(length(bv)==0) bv <- plot_pts$Outcome
    yr <- range(bv); mg <- diff(yr)*0.02; y_lim <- c(yr[1]-mg,yr[2]+mg)
  }
  p <- ggplot() +
    geom_line(data=pred_data,aes(x=Age,y=P25,color=Sex),linetype="dashed",linewidth=0.95) +
    geom_line(data=pred_data,aes(x=Age,y=P75,color=Sex),linetype="dashed",linewidth=0.95) +
    geom_line(data=pred_data,aes(x=Age,y=P50,color=Sex),linewidth=2.0,lineend="round") +
    geom_text(data=female_end,aes(x=Age+0.3,y=P25,label="P25"),color=col_f,size=2.8,hjust=0) +
    geom_text(data=female_end,aes(x=Age+0.3,y=P50,label="P50"),color=col_f,size=3.0,hjust=0,fontface="bold") +
    geom_text(data=female_end,aes(x=Age+0.3,y=P75,label="P75"),color=col_f,size=2.8,hjust=0) +
    scale_color_manual(values=c(Female=col_f,Male=col_m)) +
    scale_fill_manual( values=c(Female=col_f,Male=col_m)) +
    scale_x_continuous(expand=expansion(mult=c(0.02,0.08))) +
    scale_y_continuous(limits=y_lim,expand=c(0,0)) +
    labs(title=paste0("Normative Trajectories: ",label),
         subtitle=paste0(model_note,"\nGG distribution  |  BIC-selected fp(p,q,r)  |  Female n=",n_f,", Male n=",n_m),
         x="Age (years)",y=y_lab,color="Sex",fill="Sex",
         caption="Bold solid = P50  |  Shaded band = IQR (P25-P75)  |  Dashed = P25 / P75") +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5,size=13,face="bold"),
          plot.subtitle=element_text(hjust=0.5,size=9,color="#555",lineheight=1.4),
          plot.caption=element_text(hjust=0.5,size=8,color="#888"),
          axis.title=element_text(size=11),legend.position="top",
          panel.grid.major=element_line(color="#eeeeee"),panel.grid.minor=element_blank(),
          panel.background=element_rect(fill="white",color=NA))
  ggsave(paste0(OUTPUT_DIR,"Trajectories_",file_tag,".pdf"),plot=p,width=11,height=8,dpi=300,device=cairo_pdf)
  cat(sprintf("  ✅ Trajectories_%s.pdf\n",file_tag)); return(invisible(p))
}

# ── Fig2 三联图
plot_fig2_comprehensive <- function(pred_data,variance_data,velocity_data,plot_pts,label,file_tag,unit,y_lim=NULL,model_note="") {
  if (is.null(pred_data)||nrow(pred_data)==0){cat(sprintf("  ⚠ 跳过 Fig2 [%s]\n",file_tag));return(invisible(NULL))}
  col_f <- "#d62728"; col_m <- "#1f77b4"
  y_lab <- if(nchar(unit)>0) paste0(label," (",unit,")") else label
  if (is.null(y_lim)){av <- c(plot_pts$Outcome,pred_data$P25,pred_data$P75);yr <- range(av,na.rm=TRUE);y_lim <- c(yr[1]-0.05*diff(yr),yr[2]+0.05*diff(yr))}
  p_traj <- ggplot() +
    geom_point(data=plot_pts,aes(x=Age,y=Outcome,color=Sex),alpha=0.1,size=0.8) +
    scale_color_manual(values=c(Female=col_f,Male=col_m),guide="none") + new_scale_color() +
    geom_ribbon(data=pred_data,aes(x=Age,ymin=P25,ymax=P75,fill=Sex),alpha=0.25) +
    geom_line(data=pred_data,aes(x=Age,y=P50,color=Sex),linewidth=1.2) +
    scale_color_manual(values=c(Female=col_f,Male=col_m)) + scale_fill_manual(values=c(Female=col_f,Male=col_m)) +
    scale_y_continuous(limits=y_lim) + labs(title="Normative Trajectories",y=y_lab,x="Age (years)") +
    theme_bw()+theme(legend.position="top",plot.title=element_text(hjust=0.5,face="bold"))
  p_var <- if(!is.null(variance_data)&&nrow(variance_data)>0&&any(!is.na(variance_data$Variance)))
    ggplot(variance_data,aes(x=Age,y=Variance,color=Sex))+geom_line(linewidth=1.2)+
    scale_color_manual(values=c(Female=col_f,Male=col_m))+
    labs(title="Normative Variance",y=if(nchar(unit)>0)paste0("Variance (",unit,"\u00b2)") else "Variance",x="Age (years)")+
    theme_bw()+theme(legend.position="top",plot.title=element_text(hjust=0.5,face="bold"))
  else ggplot()+geom_blank()+labs(title="Variance (unavailable)")+theme_void()
  p_vel <- if(!is.null(velocity_data)&&nrow(velocity_data)>0&&any(!is.na(velocity_data$Velocity)))
    ggplot(velocity_data,aes(x=Age,y=Velocity,color=Sex))+geom_hline(yintercept=0,linetype="dashed",color="gray40")+
    geom_line(linewidth=1.2)+scale_color_manual(values=c(Female=col_f,Male=col_m))+
    labs(title="Rate of Growth (dy/dx)",y=if(nchar(unit)>0)paste0(unit," / year") else "unit / year",x="Age (years)")+
    theme_bw()+theme(legend.position="top",plot.title=element_text(hjust=0.5,face="bold"))
  else ggplot()+geom_blank()+labs(title="Velocity (unavailable)")+theme_void()
  p_comb <- gridExtra::grid.arrange(p_traj,p_var,p_vel,ncol=1,
                                    top=grid::textGrob(paste0("Comprehensive Analysis: ",label,"\n",model_note),gp=grid::gpar(fontface="bold",fontsize=13)))
  ggsave(paste0(OUTPUT_DIR,"Fig2_Comprehensive_",file_tag,".pdf"),plot=p_comb,width=12,height=16,dpi=300,device=cairo_pdf)
  cat(sprintf("  ✅ Fig2_Comprehensive_%s.pdf\n",file_tag)); return(invisible(p_comb))
}

# ── Fig4a 单疾病
plot_fig4a_no_stage <- function(patient_scores, label, file_tag,
                                disease_name, disease_config, output_dir,
                                return_only=FALSE) {
  df_plot <- patient_scores %>% filter(!is.na(Percentile)) %>%
    mutate(Group = factor(Sex, levels=c("Female","Male")))
  if (nrow(df_plot)==0) {
    if(!return_only) cat("  ⚠ 无有效数据，跳过 Fig4a\n")
    return(invisible(NULL))
  }
  FILLS <- c(Female=disease_config$color_female, Male=disease_config$color_male)
  sig_df <- df_plot %>% group_by(Group) %>%
    dplyr::group_modify(~ {
      pt <- .perm_test_vs_50(.x$Percentile, n_perm = 10000)
      es <- .cohens_d_vs_50(.x$Percentile)
      tibble::tibble(
        n    = sum(!is.na(.x$Percentile)),
        med  = median(.x$Percentile, na.rm = TRUE),
        pval = pt$p,
        d    = es$d,
        g    = es$g
      )
    }) %>% dplyr::ungroup() %>%
    mutate(sig   = .sig_star(pval),     # 单变量只有 2 组，FDR 没意义，用原始 p
           y_lab = 102)
  lp <- sig_df %>% select(Group,med)
  p <- ggplot(df_plot,aes(x=Group,y=Percentile,fill=Group)) +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=25,  ymax=75,  alpha=0.07,fill="#4CAF50") +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=2.5, ymax=25,  alpha=0.06,fill="#FFC107") +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=75,  ymax=97.5,alpha=0.06,fill="#FFC107") +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=0,   ymax=2.5, alpha=0.08,fill="#F44336") +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=97.5,ymax=100, alpha=0.08,fill="#F44336") +
    geom_hline(yintercept=50,linetype="dashed",linewidth=0.7,color="#555555") +
    geom_hline(yintercept=c(2.5,97.5),linetype="dotted",linewidth=0.4,color="#888888") +
    geom_hline(yintercept=c(25,75),linetype="dashed",linewidth=0.4,color="#AAAAAA") +
    geom_violin(trim=TRUE,alpha=0.55,color="white",scale="width") +
    geom_jitter(shape=21,width=0.08,size=1.4,fill="white",alpha=0.6,color="#555555",stroke=0.5) +
    geom_segment(data=lp,aes(x=Group,xend=Group,y=50,yend=med),color="white",linewidth=1.8,inherit.aes=FALSE) +
    geom_segment(data=lp,aes(x=Group,xend=Group,y=50,yend=med),color="#222222",linewidth=1.0,inherit.aes=FALSE) +
    geom_point(data=lp,aes(x=Group,y=med),shape=21,size=4.5,fill="white",color="#222222",stroke=1.5,inherit.aes=FALSE) +
    geom_text(data=lp,aes(x=Group,y=med,label=sprintf("%.0f",med)),vjust=-1.2,size=2.9,fontface="bold",color="#222222",inherit.aes=FALSE) +
    geom_text(data=sig_df,aes(x=Group,y=y_lab,label=sig),size=3.8,color="#333333",inherit.aes=FALSE) +
    geom_text(data=sig_df,aes(x=Group,y=-4,label=paste0("n=",n)),size=2.5,color="#888888",inherit.aes=FALSE) +
    scale_fill_manual(values=FILLS,guide="none") +
    scale_y_continuous(limits=c(-6,107),breaks=c(0,2.5,25,50,75,97.5,100),
                       labels=c("0","P2.5","P25","P50","P75","P97.5","100")) +
    labs(title=paste0("Centile Distributions vs Healthy Reference\n",disease_name," | ",label),
         subtitle="White line+circle=median  |  *** p<0.001  ** p<0.01  * p<0.05  ns (vs P50)",
         x=NULL,y="Centile score") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(hjust=0.5,face="bold",size=13),
          plot.subtitle=element_text(hjust=0.5,size=9,color="#555"),
          axis.text.x=element_text(size=10),
          panel.grid.major.x=element_blank(),panel.grid.minor=element_blank())
  if (!return_only) {
    ggsave(paste0(output_dir,"Fig4a_",file_tag,".pdf"),plot=p,width=8,height=7,dpi=300,device=cairo_pdf)
    cat(sprintf("  ✅ Fig4a_%s.pdf\n",file_tag))
  }
  return(invisible(p))
}

# ── Fig4c CMD 单疾病
plot_fig4c_cmd_no_stage <- function(df_all_scores, disease_name, disease_config, output_dir) {
  cmd_df <- df_all_scores %>%
    filter(!is.na(Percentile)) %>%
    group_by(PatientID, Sex, Disease) %>%
    summarise(CMD=mean(abs(Percentile-50),na.rm=TRUE), N_Var=n(), .groups="drop")
  if (nrow(cmd_df)==0) { cat("  ⚠ 无 CMD 数据\n"); return(invisible(NULL)) }
  COLS <- c(Female="#d62728", Male="#1f77b4")
  cn_ref <- 12.5
  p_main <- ggplot(cmd_df, aes(x=CMD,fill=Sex,color=Sex)) +
    geom_vline(xintercept=cn_ref, linetype="dashed", linewidth=0.8, color="#333333") +
    annotate("text", x=cn_ref+0.8, y=Inf, label="Healthy\nreference",
             hjust=0, vjust=1.3, size=3.2, color="#555555") +
    geom_density(alpha=0.40, linewidth=1.0) +
    geom_rug(alpha=0.35, length=unit(0.025,"npc")) +
    scale_fill_manual(values=COLS) +
    scale_color_manual(values=COLS) +
    scale_x_continuous(limits=c(0,50), breaks=seq(0,50,10)) +
    labs(title=paste0("CMD: Centile Mahalanobis Distance — ", disease_name),
         subtitle="CMD = mean|Percentile_j - 50| across all variables  |  larger = more atypical",
         x="CMD", y="Density") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(hjust=0.5,face="bold",size=13),
          plot.subtitle=element_text(hjust=0.5,size=9,color="#555"),
          legend.position="top", panel.grid.minor=element_blank())
  ggsave(paste0(output_dir,"Fig4c_CMD_Density.pdf"),
         plot=p_main, width=10, height=6, dpi=300, device=cairo_pdf)
  cat("  ✅ Fig4c_CMD_Density.pdf\n")
  cmd_summary <- cmd_df %>%
    group_by(Disease, Sex) %>%
    summarise(N=n(), Median_CMD=round(median(CMD,na.rm=TRUE),2),
              Mean_CMD=round(mean(CMD,na.rm=TRUE),2), SD_CMD=round(sd(CMD,na.rm=TRUE),2),
              PCT_High=round(mean(CMD>20,na.rm=TRUE)*100,1), .groups="drop")
  cat("\n--- CMD Summary (",disease_name,") ---\n",sep=""); print(cmd_summary)
  return(invisible(list(plot=p_main, data=cmd_df, summary=cmd_summary)))
}

# ★★★ v15.1：跨疾病 Fig4a — 分段布局
build_fig4a_small_panel <- function(valid_tag, var_label, disease_data_all,
                                    group_levels, fills_named) {
  df_combined <- bind_rows(lapply(names(disease_data_all), function(DIS_NAME) {
    ps <- disease_data_all[[DIS_NAME]]$all_patient_reports[[valid_tag]]
    if (is.null(ps)) return(NULL)
    ps %>% mutate(Group = paste0(DIS_NAME, "\n", Sex))
  }))
  if (is.null(df_combined) || nrow(df_combined) == 0) return(ggplot()+theme_void())
  df_combined <- df_combined %>% filter(!is.na(Percentile))
  df_combined$Group <- factor(df_combined$Group,
                              levels = group_levels[group_levels %in% unique(df_combined$Group)])
  sig_df <- df_combined %>% group_by(Group) %>%
    dplyr::group_modify(~ {
      pt <- .perm_test_vs_50(.x$Percentile, n_perm = 10000)
      es <- .cohens_d_vs_50(.x$Percentile)
      tibble::tibble(
        n    = sum(!is.na(.x$Percentile)),
        med  = median(.x$Percentile, na.rm = TRUE),
        pval = pt$p,
        d    = es$d,
        g    = es$g
      )
    }) %>% dplyr::ungroup() %>%
    mutate(sig   = .sig_star(pval),     # 每个小图只有 4 组，FDR 没意义；真正的 FDR 在 Fig5 热图里做
           y_lab = 102)
  lp <- sig_df %>% select(Group,med)
  ggplot(df_combined, aes(x=Group,y=Percentile,fill=Group)) +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=25,  ymax=75,  alpha=0.07,fill="#4CAF50") +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=2.5, ymax=25,  alpha=0.06,fill="#FFC107") +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=75,  ymax=97.5,alpha=0.06,fill="#FFC107") +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=0,   ymax=2.5, alpha=0.08,fill="#F44336") +
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=97.5,ymax=100, alpha=0.08,fill="#F44336") +
    geom_hline(yintercept=50,linetype="dashed",linewidth=0.5,color="#555555") +
    geom_hline(yintercept=c(2.5,97.5),linetype="dotted",linewidth=0.3,color="#888888") +
    geom_violin(trim=TRUE,alpha=0.55,color="white",scale="width") +
    geom_segment(data=lp,aes(x=Group,xend=Group,y=50,yend=med),
                 color="white",linewidth=1.4,inherit.aes=FALSE) +
    geom_segment(data=lp,aes(x=Group,xend=Group,y=50,yend=med),
                 color="#222222",linewidth=0.7,inherit.aes=FALSE) +
    geom_point(data=lp,aes(x=Group,y=med),shape=21,size=2.8,
               fill="white",color="#222222",stroke=1.0,inherit.aes=FALSE) +
    geom_text(data=sig_df,aes(x=Group,y=y_lab,label=sig),size=4.6,color="#333333",inherit.aes=FALSE) +
    # geom_text(data=sig_df,aes(x=Group,y=-4,label=paste0("n=",n)),size=3.2,color="#888888",inherit.aes=FALSE) +
    scale_fill_manual(values=fills_named[names(fills_named) %in% levels(df_combined$Group)],
                      guide="none") +
    scale_y_continuous(limits=c(-2,106),breaks=c(2.5,25,50,75,97.5),
                       labels=c("P2.5","P25","P50","P75","P97.5")) +
    labs(title=var_label, x=NULL, y=NULL) +
    theme_bw(base_size=14) +
    theme(plot.title=element_text(hjust=0.5,face="bold",size=12),
          axis.text.x=element_text(size=10,lineheight=1.0,angle=0),
          axis.text.y=element_text(size=10.5),
          axis.title.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin=margin(5,5,5,5))
}

build_fig4a_section_grob <- function(section, disease_data_all, label_map_full,
                                     group_levels, fills_named, ncol_panels = 6) {
  use_tags <- get_section_selected_tags(section)
  common_tags <- Reduce(intersect, lapply(disease_data_all,
                                          function(d) names(d$all_patient_reports)))
  use_tags <- use_tags[use_tags %in% common_tags]
  if (length(use_tags) == 0) {
    body <- grid::textGrob(sprintf("(%s) 该节无可用变量", section$label),
                           gp=grid::gpar(fontsize=9, col="#888"))
    header <- grid::textGrob(section$label, x=unit(0.01,"npc"), hjust=0,
                             gp=grid::gpar(fontface="bold", fontsize=12, col="#1B3A5C"))
    divider <- grid::linesGrob(x=unit(c(0.005,0.995),"npc"), y=unit(c(0.5,0.5),"npc"),
                               gp=grid::gpar(col="#1B3A5C", lwd=1.2))
    return(list(header=header, divider=divider, body=body, nrow=1))
  }
  panels <- lapply(use_tags, function(tag) {
    var_label <- if (!is.na(label_map_full[tag])) label_map_full[tag] else tag
    build_fig4a_small_panel(tag, var_label, disease_data_all, group_levels, fills_named)
  })
  n      <- length(panels)
  nrow_p <- ceiling(n / ncol_panels)
  if (n < nrow_p * ncol_panels) {
    for (k in (n+1):(nrow_p*ncol_panels)) panels[[k]] <- ggplot()+theme_void()
  }
  panel_grid <- gridExtra::arrangeGrob(grobs=panels, ncol=ncol_panels, nrow=nrow_p)
  header <- grid::textGrob(section$label, x=unit(0.01,"npc"), hjust=0,
                           gp=grid::gpar(fontface="bold", fontsize=12, col="#1B3A5C"))
  divider <- grid::linesGrob(x=unit(c(0.005,0.995),"npc"), y=unit(c(0.5,0.5),"npc"),
                             gp=grid::gpar(col="#1B3A5C", lwd=1.2))
  list(header=header, divider=divider, body=panel_grid, nrow=nrow_p)
}

plot_fig4a_sectioned_multi_disease <- function(disease_data_all, SECTIONS,
                                               output_dir,
                                               ncol_panels = SELECTED_SUMMARY_NCOL,
                                               row_height = 2.0,
                                               header_height = 0.18,
                                               divider_height = 0.02) {
  if (length(disease_data_all) == 0) {
    cat("  ⚠ 无疾病数据，跳过 Fig4a_Selected_Summary_MultiDisease\n"); return(invisible(NULL))
  }
  
  label_map_full <- setNames(sapply(TARGETS, `[[`, "label"),
                             sapply(TARGETS, `[[`, "file_tag"))
  
  COL_FEMALE <- "#d62728"; COL_MALE <- "#1f77b4"
  group_levels <- character(0); fills_named <- character(0)
  for (DIS_NAME in names(disease_data_all)) {
    fkey <- paste0(DIS_NAME,"\nFemale"); mkey <- paste0(DIS_NAME,"\nMale")
    fills_named[fkey] <- COL_FEMALE; fills_named[mkey] <- COL_MALE
    group_levels <- c(group_levels, fkey, mkey)
  }
  
  section_grobs <- lapply(SECTIONS, build_fig4a_section_grob,
                          disease_data_all = disease_data_all,
                          label_map_full   = label_map_full,
                          group_levels     = group_levels,
                          fills_named      = fills_named,
                          ncol_panels      = ncol_panels)
  
  final_list  <- list(); heights_vec <- numeric(0)
  for (sg in section_grobs) {
    final_list  <- c(final_list, list(sg$header, sg$divider, sg$body))
    heights_vec <- c(heights_vec, header_height, divider_height, sg$nrow * row_height)
  }
  
  total_vars <- sum(sapply(SECTIONS, function(s) length(get_section_selected_tags(s))))
  pdf_path <- paste0(output_dir, "Fig4a_Selected_Summary_MultiDisease.pdf")
  disease_str <- paste(names(disease_data_all), collapse=" vs ")
  p_all <- gridExtra::grid.arrange(
    grobs   = final_list,
    ncol    = 1,
    heights = heights_vec,
    bottom  = grid::textGrob(
      "White line+circle = group median  |  Dashed = P50 healthy reference  |  Female = red, Male = blue  |  *** p<0.001  ** p<0.01  * p<0.05  ns (vs P50)",
      gp = grid::gpar(fontsize=15, col="#555"))
  )
  total_height <- sum(heights_vec) + 1.2
  save_key_figure(p_all, gsub("\\.pdf$","",pdf_path), width=16, height=total_height)
}

# ★★★ v15.1：跨疾病 显著性热图（替代旧的相关性热图）
# ★★★ v15.2 通用版：跨疾病显著性热图（矩形 + 环形 双输出）
#   part_label：图正中央显示的文字（"Complexity" / "Tortuosity" / "Caliber"）
#   ring_mode： TRUE  = 同时输出矩形版 + 环形版（推荐 Complexity / Caliber）
#               FALSE = 只输出矩形版（Tortuosity 参数太少，环形不好看）
plot_significance_heatmap_multi_disease <- function(disease_data_all, SECTIONS, output_dir,
                                                    part_label = "",
                                                    ring_mode = TRUE) {
  if (length(disease_data_all) == 0) {
    cat("  ⚠ 无疾病数据，跳过 显著性热图\n"); return(invisible(NULL))
  }
  
  label_map_full <- setNames(sapply(TARGETS, `[[`, "label"),
                             sapply(TARGETS, `[[`, "file_tag"))
  
  # ───────── 1. 构建长表 ─────────
  rows <- list()
  for (sec_i in seq_along(SECTIONS)) {
    sec <- SECTIONS[[sec_i]]
    use_tags <- get_section_selected_tags(sec)
    for (tag in use_tags) {
      var_lab <- if (!is.na(label_map_full[tag])) label_map_full[tag] else tag
      for (DIS_NAME in names(disease_data_all)) {
        ps <- disease_data_all[[DIS_NAME]]$all_patient_reports[[tag]]
        if (is.null(ps)) next
        for (sx in c("Female","Male")) {
          v <- ps$Percentile[ps$Sex == sx]
          v <- v[!is.na(v)]
          n_v <- length(v)
          if (n_v == 0) {
            rows[[length(rows)+1]] <- data.frame(
              Section_Idx=sec_i, Section=sec$label, Tag=tag, Variable=var_lab,
              Disease=DIS_NAME, Sex=sx, N=0,
              Median_Pct=NA_real_, Median_Diff=NA_real_,
              Mean_Pct=NA_real_,   Mean_Diff=NA_real_,
              Cohens_d=NA_real_,   Hedges_g=NA_real_,
              Pvalue=NA_real_,     Pvalue_FDR=NA_real_, Sig="", stringsAsFactors=FALSE)
            next
          }
          med      <- median(v); med_diff <- med - 50
          pt       <- if (n_v >= 5) .perm_test_vs_50(v, n_perm = 10000)
          else list(obs = NA_real_, p = NA_real_)
          es       <- if (n_v >= 5) .cohens_d_vs_50(v)
          else list(d = NA_real_, g = NA_real_)
          rows[[length(rows)+1]] <- data.frame(
            Section_Idx=sec_i, Section=sec$label, Tag=tag, Variable=var_lab,
            Disease=DIS_NAME, Sex=sx, N=n_v,
            Median_Pct=round(med,2), Median_Diff=round(med_diff,2),
            Mean_Pct=round(mean(v),2), Mean_Diff=round(pt$obs,2),
            Cohens_d=round(es$d,3), Hedges_g=round(es$g,3),
            Pvalue=signif(pt$p,3), Pvalue_FDR=NA_real_, Sig="",   # FDR/Sig 在 long_df 全表后统一算
            stringsAsFactors=FALSE)
        }
      }
    }
  }
  if (length(rows) == 0) {
    cat("  ⚠ 无可用数据，跳过 显著性热图\n"); return(invisible(NULL))
  }
  long_df <- do.call(rbind, rows)
  # ★ Nature SI 10.1：跨所有 disease × sex × variable 做全局 BH FDR
  long_df$Pvalue_FDR <- p.adjust(long_df$Pvalue, method = "BH")
  long_df$Sig        <- .sig_star(long_df$Pvalue_FDR)
  
  long_df$DiseaseSex <- paste0(long_df$Disease, "\n", long_df$Sex)
  ds_levels <- unlist(lapply(names(disease_data_all),
                             function(d) paste0(d, "\n", c("Female","Male"))))
  long_df$DiseaseSex <- factor(long_df$DiseaseSex, levels = ds_levels)
  
  # ★ 变量顺序：严格按 SECTIONS 顺序排列
  ordered_tags <- unlist(lapply(SECTIONS, function(s) get_section_selected_tags(s)))
  ordered_tags <- ordered_tags[ordered_tags %in% long_df$Tag]
  ordered_vars <- unique(label_map_full[ordered_tags])
  ordered_vars <- ordered_vars[!is.na(ordered_vars)]
  
  long_df$cell_label <- ifelse(is.na(long_df$Median_Pct), "",
                               paste0(sprintf("%.0f", long_df$Median_Pct),
                                      ifelse(nchar(long_df$Sig)>0 & long_df$Sig!="ns",
                                             paste0("\n", long_df$Sig), "")))
  
  fill_lim <- c(-50, 50)
  
  # ───────── 2. 矩形版热图（保留原版） ─────────
  long_df_rect <- long_df
  long_df_rect$Variable <- factor(long_df_rect$Variable, levels = rev(ordered_vars))
  
  p_heat_rect <- ggplot(long_df_rect, aes(x = DiseaseSex, y = Variable, fill = Median_Diff)) +
    geom_tile(color="white", linewidth=0.4) +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                         midpoint=0, limits=fill_lim, oob=scales::squish,
                         name="median(Percentile)\n      - 50",
                         na.value="#EEEEEE") +
    geom_text(aes(label=cell_label), size=3.6, color="black", lineheight=0.85) +
    labs(title=sprintf("Per-variable Centile Deviation from Healthy P50 — %s",
                       gsub("\n", " ", part_label)),
         subtitle="Color = median(Percentile) - 50  |  Red = HIGH, Blue = LOW  |  *** p<0.001  ** p<0.01  * p<0.05",
         x=NULL, y=NULL) +
    theme_minimal(base_size=14) +
    theme(plot.title    = element_text(hjust=0.5, face="bold", size=18),
          plot.subtitle = element_text(hjust=0.5, size=12, color="#555"),
          axis.text.x   = element_text(size=14, lineheight=0.9),
          axis.text.y   = element_text(size=11),
          panel.grid    = element_blank(),
          legend.position = "right",
          legend.title  = element_text(size=12),
          legend.text   = element_text(size=11))
  
  pdf_rect <- paste0(output_dir, "Fig5_Significance_Heatmap_MultiDisease.pdf")
  n_var <- length(ordered_vars); n_col <- length(ds_levels)
  width_rect  <- max(8, 2.5 + 0.95 * n_col + 2.5)
  height_rect <- max(6, 2.0 + 0.22 * n_var)
  save_key_figure(p_heat_rect, gsub("\\.pdf$","",pdf_rect), width=width_rect, height=height_rect)
  cat(sprintf("  ✅ Fig5_Significance_Heatmap_MultiDisease.pdf (矩形, %d 行 × %d 列)\n",
              n_var, n_col))
  
  # ───────── 3. 环形版热图（参数 ≥ 6 才出，否则跳过） ─────────
  if (ring_mode && n_var >= 6) {
    
    # ★ 节间插入空"哨兵"扇形作为视觉分隔，并记录哨兵的 x 位置（用于放置 ring 标签）
    sentinel_vars <- character(0)
    long_df_ring  <- long_df
    long_df_ring$Variable <- as.character(long_df_ring$Variable)
    
    augmented_vars       <- character(0)
    augmented_sections   <- integer(0)
    sentinel_x_positions <- integer(0)
    for (sec_i in seq_along(SECTIONS)) {
      sec <- SECTIONS[[sec_i]]
      sec_tags <- get_section_selected_tags(sec)
      sec_tags <- sec_tags[sec_tags %in% long_df$Tag]
      sec_vars <- unname(label_map_full[sec_tags])
      sec_vars <- sec_vars[!is.na(sec_vars)]
      
      augmented_vars     <- c(augmented_vars, sec_vars)
      augmented_sections <- c(augmented_sections, rep(sec_i, length(sec_vars)))
      
      # 每个模块后都加一个空白扇区：
      # 前几个是模块间缺口，最后一个是首尾闭合处缺口
      sentinel <- paste0("__sep_", sec_i, "__")
      sentinel_vars <- c(sentinel_vars, sentinel)
      augmented_vars     <- c(augmented_vars, sentinel)
      augmented_sections <- c(augmented_sections, NA_integer_)
      sentinel_x_positions <- c(sentinel_x_positions, length(augmented_vars))
    }
    
    sentinel_rows <- expand.grid(Variable = sentinel_vars,
                                 DiseaseSex = ds_levels,
                                 stringsAsFactors = FALSE)
    sentinel_rows$Median_Diff <- NA_real_
    sentinel_rows$cell_label  <- ""
    sentinel_rows$Sig         <- ""
    long_df_ring <- bind_rows(long_df_ring, sentinel_rows)
    
    long_df_ring$Variable   <- factor(long_df_ring$Variable, levels = augmented_vars)
    long_df_ring$DiseaseSex <- factor(long_df_ring$DiseaseSex, levels = ds_levels)
    long_df_ring$RingLevel  <- as.integer(long_df_ring$DiseaseSex)
    
    n_rings  <- length(ds_levels)
    inner_r  <- 1.5      # 中央空白半径
    n_total  <- length(augmented_vars)
    
    # ★★★ 修复1：变量标签沿半径外读（避免重叠）
    real_idx <- which(!grepl("^__sep_", augmented_vars))
    var_label_df <- data.frame(
      x_pos      = real_idx,
      label_text = augmented_vars[real_idx],
      section_id = augmented_sections[real_idx],
      stringsAsFactors = FALSE
    )
    # 数学角度：12 点 = 90°；保留原来的 90 - ...，只修正 12 点边界
    .raw_va <- 90 + (var_label_df$x_pos - 1) * 360 / n_total
    var_label_df$math_angle <- ((.raw_va + 180) %% 360) - 180
    
    # 左半圈翻转使字正向阅读
    # 关键修改：把正上方 90° 也归入翻转侧，避免 12 点方向朝向不统一
    # 先按常规规则判断
    var_label_df$flip <- with(var_label_df, math_angle > 90 | math_angle < -90)
    
    # 12 点附近的标签不要单独判断，跟随同一模块最近的非 12 点标签
    top_eps <- 0.5
    top_idx <- which(abs(var_label_df$math_angle - 90) < top_eps)
    
    for (i in top_idx) {
      same_sec <- which(
        var_label_df$section_id == var_label_df$section_id[i] &
          abs(var_label_df$math_angle - 90) >= top_eps
      )
      
      if (length(same_sec) > 0) {
        nearest <- same_sec[which.min(abs(var_label_df$x_pos[same_sec] - var_label_df$x_pos[i]))]
        var_label_df$flip[i] <- var_label_df$flip[nearest]
      }
    }
    
    var_label_df$angle <- with(var_label_df, ifelse(flip, math_angle + 180, math_angle))
    var_label_df$angle <- ((var_label_df$angle + 180) %% 360) - 180
    var_label_df$hjust <- with(var_label_df, ifelse(flip, 1, 0))
    
    # ★★★ 修复2：Disease×Sex 标签放在第一个哨兵位置（节间空白处），沿半径方向旋转
    if (length(sentinel_x_positions) > 0) {
      ring_x_pos <- sentinel_x_positions[1]
    } else {
      ring_x_pos <- 1   # 退化情况（理论上 ring_mode 已要求 n_var >= 6 + 多节）
    }
    .raw_ring <- 90 + (ring_x_pos - 1) * 360 / n_total
    ring_math_angle <- ((.raw_ring + 180) %% 360) - 180
    
    ring_flip       <- ring_math_angle >= 90 | ring_math_angle < -90
    ring_angle      <- ifelse(ring_flip, ring_math_angle + 180, ring_math_angle)
    ring_angle      <- ((ring_angle + 180) %% 360) - 180
    
    ring_label_df <- data.frame(
      x_pos      = ring_x_pos,
      RingLevel  = seq_along(ds_levels),
      label_text = ds_levels,
      angle      = ring_angle,
      stringsAsFactors = FALSE
    )
    
    # ★★★ 修复3：移除外圈节标题（用户不需要）
    
    p_heat_ring <- ggplot(long_df_ring,
                          aes(x = Variable, y = RingLevel + inner_r, fill = Median_Diff)) +
      geom_tile(color="white", linewidth=0.3) +
      scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                           midpoint=0, limits=fill_lim, oob=scales::squish,
                           name="median(Pct)\n   - 50",
                           na.value="#FFFFFF") +
      geom_text(aes(label=cell_label), size=3.4, color="black", lineheight=0.85) +
      # ★ 变量标签：沿半径外读
      geom_text(data = var_label_df,
                aes(x = x_pos, y = n_rings + inner_r + 0.4,
                    label = label_text,
                    angle = angle, hjust = hjust),
                inherit.aes = FALSE,
                size = 4.0, color = "#222222") +
      # ★ Disease × Sex 标签：嵌在节间哨兵空白处，沿半径方向，字号加大
      geom_text(data = ring_label_df,
                aes(x = x_pos, y = RingLevel + inner_r,
                    label = label_text,
                    angle = angle),
                inherit.aes = FALSE,
                hjust = 0.5, vjust = 0.5,
                size = 3.6, color = "#1B3A5C", fontface = "bold",
                lineheight = 0.85) +
      # 中央 part 标题
      annotate("text", x = 0.5, y = 0,
               label = part_label,
               size = 9, fontface = "bold", color = "#1B3A5C") +
      scale_y_continuous(limits = c(0, n_rings + inner_r + 2.5)) +
      coord_polar(theta = "x", start = -pi/length(augmented_vars), direction = -1, clip = "off") +
      labs(title    = sprintf("Per-variable Centile Deviation from Healthy P50 — %s (Circular)",
                              part_label),
           subtitle = "Sector = variable  |  Ring (inner→outer) = Disease × Sex  |  Red = HIGH, Blue = LOW",
           x=NULL, y=NULL) +
      theme_minimal(base_size=13) +
      theme(plot.title    = element_text(hjust=0.5, face="bold", size=19),
            plot.subtitle = element_text(hjust=0.5, size=12, color="#555"),
            axis.text.x   = element_blank(),     # ← 关闭默认变量标签（用 geom_text 替代）
            axis.text.y   = element_blank(),     # ← 关闭默认 ring 标签（用 geom_text 替代）
            axis.ticks    = element_blank(),
            panel.grid    = element_blank(),
            legend.position = "right",
            plot.margin   = margin(40, 40, 40, 40))
    
    pdf_ring <- paste0(output_dir, "Fig5_Significance_Heatmap_Circular.pdf")
    side <- max(13, min(20, 10 + n_var * 0.12))  # 边距更大，给外圈标签留空间
    save_key_figure(p_heat_ring, gsub("\\.pdf$","",pdf_ring), width=side, height=side)
    cat(sprintf("  ✅ Fig5_Significance_Heatmap_Circular.pdf (环形, %d sectors × %d rings)\n",
                n_var, n_rings))
    
  } else if (ring_mode && n_var < 6) {
    cat(sprintf("  ⓘ 参数数量 %d < 6，跳过环形版（参数太少环形图不好看）\n", n_var))
  }
  
  # ───────── 4. Excel 输出（不变） ─────────
  long_out <- long_df %>%
    select(Section, Variable, Tag, Disease, Sex, N,
           Median_Pct, Median_Diff, Mean_Pct, Mean_Diff,
           Cohens_d, Hedges_g, Pvalue, Pvalue_FDR, Sig)
  wide_med  <- long_df %>%
    mutate(DSstr = paste(Disease, Sex, sep="_")) %>%
    select(Section, Variable, DSstr, Median_Diff) %>%
    tidyr::pivot_wider(names_from=DSstr, values_from=Median_Diff)
  wide_d    <- long_df %>%                                  # ★ 新增
    mutate(DSstr = paste(Disease, Sex, sep="_")) %>%
    select(Section, Variable, DSstr, Cohens_d) %>%
    tidyr::pivot_wider(names_from=DSstr, values_from=Cohens_d)
  wide_pval <- long_df %>%
    mutate(DSstr = paste(Disease, Sex, sep="_")) %>%
    select(Section, Variable, DSstr, Pvalue_FDR) %>%        # ★ 改为 FDR 后的 p
    tidyr::pivot_wider(names_from=DSstr, values_from=Pvalue_FDR)
  wide_sig  <- long_df %>%
    mutate(DSstr = paste(Disease, Sex, sep="_")) %>%
    select(Section, Variable, DSstr, Sig) %>%
    tidyr::pivot_wider(names_from=DSstr, values_from=Sig)
  
  wb_sig <- createWorkbook()
  hs_sig <- createStyle(fontColour="#FFFFFF", bgFill="#5D4E60", halign="center",
                        textDecoration="bold", border="Bottom")
  add_sh <- function(wb,nm,dat){
    addWorksheet(wb,nm); writeData(wb,nm,as.data.frame(dat))
    addStyle(wb,nm,hs_sig,rows=1,cols=1:ncol(dat),gridExpand=TRUE)
    setColWidths(wb,nm,cols=1:ncol(dat),widths="auto")
  }
  add_sh(wb_sig, "Long_All",        long_out)
  add_sh(wb_sig, "Wide_MedianDiff", wide_med)
  add_sh(wb_sig, "Wide_CohensD",    wide_d)        # ★ 新增
  add_sh(wb_sig, "Wide_Pvalue",     wide_pval)
  add_sh(wb_sig, "Wide_Sig",        wide_sig)
  saveWorkbook(wb_sig, paste0(output_dir, "Significance_Matrix_MultiDisease.xlsx"),
               overwrite=TRUE)
  cat("  ✅ Significance_Matrix_MultiDisease.xlsx (Long_All / Wide_MedianDiff / Wide_Pvalue / Wide_Sig)\n")
  
  return(invisible(long_df))
}

# ── 汇总小图（健康轨迹）— 强制无图例
make_small_traj_panel <- function(result) {
  if (is.null(result)) return(ggplot()+theme_void())
  pred_data <- result$pred_data; plot_pts <- result$plot_pts
  label <- result$label; unit <- result$unit
  if (is.null(pred_data)||nrow(pred_data)==0) return(ggplot()+theme_void())
  col_f <- "#d62728"; col_m <- "#1f77b4"
  title_lab <- if(nchar(unit)>0) paste0(label," (",unit,")") else label
  bv <- c(pred_data$P25,pred_data$P75); bv <- bv[!is.na(bv)]
  if(length(bv)==0) bv <- plot_pts$Outcome
  yr <- range(bv); mg <- diff(yr)*0.03
  y_lim <- c(yr[1]-mg,yr[2]+mg)
  ggplot() +
    geom_line(data=pred_data,aes(x=Age,y=P25,color=Sex),linetype="dashed",linewidth=0.95) +
    geom_line(data=pred_data,aes(x=Age,y=P75,color=Sex),linetype="dashed",linewidth=0.95) +
    geom_line(data=pred_data,aes(x=Age,y=P50,color=Sex),linewidth=2.0,lineend="round") +
    scale_color_manual(values=c(Female="#d62728",Male="#1f77b4")) +
    scale_y_continuous(limits=y_lim,expand=c(0,0)) +
    scale_x_continuous(breaks=c(45,55,65),labels=c("45","55","65yr"),
                       expand=expansion(mult=c(0.02,0.08))) +
    coord_cartesian(xlim=c(NA, 67), clip="off") +
    labs(title=title_lab, x=NULL, y=NULL) +
    theme_bw(base_size=14) +
    theme(plot.title=element_text(hjust=0.5,size=12,face="bold",margin=margin(b=3)),
          axis.title=element_blank(),
          axis.text=element_text(size=11),
          legend.position="none",
          panel.grid.major=element_line(color="#f0f0f0",linewidth=0.4),
          panel.grid.minor=element_blank(),plot.margin=margin(6,6,6,6))
}

make_legend_grob <- function() {
  col_f <- "#d62728"; col_m <- "#1f77b4"
  dummy <- data.frame(Age=c(40,40,60,60), P50=c(1,2,1,2),
                      Sex=factor(c("Female","Male","Female","Male"),levels=c("Female","Male")))
  p_dummy <- ggplot(dummy, aes(x=Age,y=P50,color=Sex,fill=Sex)) +
    geom_line(linewidth=2.0) +
    scale_color_manual(values=c(Female=col_f,Male=col_m)) +
    scale_fill_manual( values=c(Female=col_f,Male=col_m)) +
    theme_bw(base_size=10) +
    theme(legend.position="bottom",legend.title=element_text(face="bold"),legend.text=element_text(size=10))
  cowplot::get_legend(p_dummy)
}

# 构建带小标题的分段大图（每节 6 列；不足列自动填空白）
build_section_grob <- function(section, results_list, ncol_panels = 6) {
  use_tags <- get_section_selected_tags(section)                          # ★ 新增
  panels   <- lapply(use_tags, function(tag) make_small_traj_panel(results_list[[tag]]))  # ★ 改用 use_tags
  n        <- length(panels)
  nrow_p   <- ceiling(n / ncol_panels)
  if (n < nrow_p * ncol_panels) {
    for (k in (n+1):(nrow_p*ncol_panels)) panels[[k]] <- ggplot()+theme_void()
  }
  panel_grid <- gridExtra::arrangeGrob(grobs=panels, ncol=ncol_panels, nrow=nrow_p)
  
  header <- grid::textGrob(
    section$label,
    x = unit(0.01, "npc"), hjust = 0,
    gp = grid::gpar(fontface="bold", fontsize=12, col="#1B3A5C")
  )
  divider <- grid::linesGrob(
    x = unit(c(0.005, 0.995), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    gp = grid::gpar(col="#1B3A5C", lwd=1.2)
  )
  
  list(header=header, divider=divider, body=panel_grid, nrow=nrow_p)
}

build_full_layout_with_sections <- function(SECTIONS, results_list,
                                            ncol_panels=6,
                                            row_height=2.6,
                                            header_height=0.55,
                                            divider_height=0.06) {
  section_grobs <- lapply(SECTIONS, build_section_grob, results_list, ncol_panels)
  
  final_list  <- list()
  heights_vec <- numeric(0)
  for (sg in section_grobs) {
    final_list  <- c(final_list, list(sg$header, sg$divider, sg$body))
    heights_vec <- c(heights_vec, header_height, divider_height, sg$nrow * row_height)
  }
  list(grobs=final_list, heights=heights_vec, total_nrow=sum(sapply(section_grobs,`[[`,"nrow")))
}

# ========================================================================================
# 6. 画图阶段
# ========================================================================================
cat("\n>>> 阶段3: 生成所有 PDF 图 <<<\n")

# 6.1 健康参考图
for (i in seq_along(TARGETS)) {
  tgt <- TARGETS[[i]]; res <- ALL_RESULTS[[tgt$file_tag]]
  if (is.null(res)) next
  cat(sprintf("\n  === [%d/%d] %s 健康图 ===\n",i,length(TARGETS),tgt$label))
  if (!is.null(res$df_h)) plot_violins(res$df_h, tgt$label, tgt$file_tag, tgt$unit)
  plot_trajectories(res$pred_data, res$plot_pts, tgt$label, tgt$file_tag, tgt$unit,
                    NULL, res$model_note)
  plot_fig2_comprehensive(res$pred_data, res$variance_data, res$velocity_data,
                          res$plot_pts, tgt$label, tgt$file_tag, tgt$unit,
                          NULL, res$model_note)
}

# 6.2 各疾病单独 Fig4a / Fig4c
catch_plot_warnings <- function(expr, context_label) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("fewer than two datapoints|drop\\s*=\\s*FALSE", msg, ignore.case = TRUE)) {
        cat(sprintf("    ⚠ [%s] 触发警告: %s\n", context_label, trimws(msg)))
        invokeRestart("muffleWarning")
      }
    }
  )
}

if (length(disease_data_all) > 0) {
  cat("\n--- 各疾病单变量 Fig4a / Fig4c ---\n")
  for (DIS_NAME in names(disease_data_all)) {
    DIS_DATA   <- disease_data_all[[DIS_NAME]]
    DISEASE_DIR <- paste0(OUTPUT_DIR, DIS_NAME, "/")
    if (!dir.exists(DISEASE_DIR)) dir.create(DISEASE_DIR, recursive=TRUE)
    cat(sprintf("\n  >>> 疾病: %s\n", DIS_NAME))
    
    label_map_full <- setNames(sapply(TARGETS, `[[`, "label"),
                               sapply(TARGETS, `[[`, "file_tag"))
    for (tag in names(DIS_DATA$all_patient_reports)) {
      var_label <- label_map_full[tag]
      if (is.na(var_label)) var_label <- tag
      catch_plot_warnings(
        plot_fig4a_no_stage(DIS_DATA$all_patient_reports[[tag]], var_label, tag,
                            DIS_NAME, DIS_DATA$config, DISEASE_DIR),
        context_label = sprintf("Fig4a / %s / %s", DIS_NAME, tag)
      )
    }
    
    if (!is.null(DIS_DATA$df_all_scores)) {
      cat(sprintf("\n  --- Fig4c CMD (%s) ---\n", DIS_NAME))
      catch_plot_warnings(
        plot_fig4c_cmd_no_stage(DIS_DATA$df_all_scores, DIS_NAME, DIS_DATA$config, DISEASE_DIR),
        context_label = sprintf("Fig4c_CMD_Density / %s", DIS_NAME)
      )
    }
  }
  
  # 6.3 跨疾病合并 Fig4a — 分段布局（与 Trajectories 大图一致；selected 由 SECTIONS 控制）
  cat("\n--- 跨疾病 Fig4a_Selected_Summary（分段, 全参数; CKD vs AD）---\n")
  catch_plot_warnings(
    plot_fig4a_sectioned_multi_disease(disease_data_all, SECTIONS, OUTPUT_DIR),
    context_label = "Fig4a_Selected_Summary_MultiDisease"
  )
  
  # 6.4 ★★★ v15.1 重写：跨疾病显著性热图
  cat("\n--- 跨疾病 Fig5_Significance_Heatmap_MultiDisease ---\n")
  catch_plot_warnings(
    plot_significance_heatmap_multi_disease(disease_data_all, SECTIONS, OUTPUT_DIR,
                                            part_label = "Caliber",
                                            ring_mode  = TRUE),
    context_label = "Fig5_Significance_Heatmap_MultiDisease"
  )
}

# 6.5 ★★★ v15.0：分段大图（4 节，每节 6 列，41 参数）
cat("\n--- 汇总 Trajectories 大图（4 节, 每节 6 列, 41 参数）---\n")

layout_obj <- build_full_layout_with_sections(SECTIONS, ALL_RESULTS,
                                              ncol_panels = 6,
                                              row_height = 1.8,
                                              header_height = 0.18,
                                              divider_height = 0.02)

shared_legend <- make_legend_grob()
final_grobs   <- c(layout_obj$grobs, list(shared_legend))
final_heights <- c(layout_obj$heights, 0.45)

pdf_all <- paste0(OUTPUT_DIR, "Part3_Caliber_Trajectories_Sectioned.pdf")
p_all <- gridExtra::grid.arrange(
  grobs   = final_grobs,
  ncol    = 1,
  heights = final_heights,
  bottom  = grid::textGrob(
    "Bold solid = P50  |  Dashed = P25 / P75  |  Red = Female, Blue = Male",
    gp = grid::gpar(fontsize=12, col="#555"))
)
total_height <- sum(final_heights) + 1.2
save_key_figure(p_all, gsub("\\.pdf$","",pdf_all), width=15, height=total_height)
cat(sprintf("  ✅ %s\n", basename(pdf_all)))

cat("\n",strrep("=",80),"\n")
cat("  ✅ Part 3 Caliber v15.1 全部完成！\n")
cat("  输出目录:", OUTPUT_DIR, "\n")
cat("  ★ v15.1 新功能：\n")
cat("    - 4 节分段大图（Trajectories）：每节 6 列、中文小标题\n")
cat("    - 简洁短标签（Diameter Artery_L1 / Diameter Artery_L1_Top1 / AVR_L1 等）\n")
cat("    - Fig4a_Selected_Summary_MultiDisease 重写为分段布局（默认全 41 参数）\n")
cat("    - Fig5_Significance_Heatmap_MultiDisease：分性别 × 疾病 的偏离量与显著性热图\n")
cat("    - Significance_Matrix_MultiDisease.xlsx (Long_All / Wide_MedianDiff / Wide_Pvalue / Wide_Sig)\n")
cat("    - 保留 ensure_avr()：自动从 top6 直径计算 AVR_1/2/3\n")
cat(strrep("=",80),"\n\n")