snp_scan36_obj <- scan1(snpprob_collapse(snpprobs36_obj,
                                         snp_action), 
                        phe_df, K_chr, cov_mx)

top_snps36_tbl <- get_top_snps_tbl(snp_scan36_obj)

(patterns <- summary(top_snps36_tbl) %>%
    filter(max_lod >= 3) %>%
    mutate(contrast = snp_action) %>%
    mutate(max_snp = "") %>%
    arrange(desc(max_lod)))

scan_pat <- 
  scan_pattern(probs36_obj,
               phe_df[,phe_par$pheno_names, drop=FALSE],
               K_chr, cov_mx, patterns,
               haplos, diplos)

#(pattern <- sdp_to_pattern(patterns$sdp))

(pattern_cont <- 
    (scan_pat$patterns %>%
       filter(sdp_to_pattern(sdp) == pattern))$contrast)

## Produce plot
pdf(paste0("images/snp_allele_3_",
           pattern, "_",
           snp_action, ".pdf"),
    width=10,height=6)
par(mfrow=c(1,2))
plot(scan_pat, "coef", pattern_cont)
plot.new()
vps <- baseViewports()
pushViewport(vps$figure)
vp1 <-plotViewport(c(1,1,1,1)) 
print(plot(scan_pat, "lod", pattern_cont),
      vp = vp1)
dev.off()

