library(data.table)


modules <- fread("/home/ttamaki/New/results/modules/exp_mod_r7.6_c0.16.tsv")

conj <- fread("/home/ttamaki/New/results/CONJScan/read_best_solution.tsv")

t4cp_genes <- conj[gene_name %in% c("T4SS_t4cp1","T4SS_t4cp2","T4SS_tcpA")]$hit_id

setnames(modules, old = "gene", new = "gene_n")
modules

modules[ , gene := paste0(replicon,"_",gene_n)]
modules
t4cp_communities <- unique(modules[gene %in% t4cp_genes$hit_id]$community_id)
t4cp_modules <- modules[community_id %in% t4cp_communities]
annot <- fread("/home/ttamaki/New/results/member_cluster_id_annotated.tsv")
unique(annot[cluster_n %in% t4cp_communities]$signature_description)

fwrite(t4cp_modules, "/home/ttamaki/New/results/modules/t4cp_modules_r7.6_c0.16.tsv", sep = "\t")



### T4cp fisher
library(data.table)
library(tidyr)

# Paths
CONJScan_path <- "/home/ttamaki/New/results/CONJScan/systems_dt_T4SS.tsv"
member_id_path <- "/home/ttamaki/New/results/cooccurrence_out_member_id.tsv"
clusters_communities_path <- "/home/ttamaki/New/results/communities/exp_com_r7.6_c0.16.tsv" 

# Clusters
member_id <- fread(member_id_path, header = TRUE)

# Comunidades
clusters_communities <- fread(clusters_communities_path, header = TRUE)
clusters_communities[,gene := paste0(id,"_",gene_n)]
clusters_communities

clusters_communities <- merge(
  clusters_communities,
  member_id,
  by = "gene",
)

communities <- merge(
    member_id,
    clusters_communities[, .(cluster_n, community_n)],
    by = "cluster_n",
    all.x = TRUE
)

communities <- na.omit(communities)
setnames(communities, old = c("gene", "cluster_n", "community_n"), new = c("gene", "cluster_id", "community_id"))
communities[,c("plasmid","gene_n"):= tstrsplit(gene, "_", fixed=TRUE)]
communities[, t4cp := any(gene %in% t4cp_genes), by = community_id]
t4cp_communities <- communities[t4cp == TRUE]
t4cp_communities[, t4cp := NULL]

### Adiciona o id do módulo
t4cp_communities[
    ,
    module_id := .GRP,
    by = .(community_id, plasmid)
]

### Tira módulos que possuem apenas um gene
t4cp_communities <- t4cp_communities[
    t4cp_communities[, .I[.N > 1], by = module_id]$V1
]

t4cp_communities[, in_t4cp := gene %in% t4cp_genes]

module_info <- t4cp_communities[
    ,
    .(t4cp = any(in_t4cp)),
    by = .(community_id, module_id)
]


module_cluster <- unique(
    t4cp_communities[
        ,
        .(community_id, module_id, cluster_id)
    ]
)

module_cluster <- merge(
    module_cluster,
    module_info,
    by = c("community_id", "module_id")
)


clusters <- unique(module_cluster$cluster_id)

resultado <- rbindlist(lapply(clusters, function(cl){

    presente <- unique(
        module_cluster[
            ,
            .(
                community_id,
                module_id,
                t4cp,
                presente = cluster_id == cl
            )
        ]
    )

    presente <- presente[
        ,
        .(presente = any(presente)),
        by = .(community_id, module_id, t4cp)
    ]

    ft <- fisher.test(table(presente$t4cp, presente$presente),
                      alternative = "greater")

    data.table(
        cluster_id = cl,
        pvalue = ft$p.value,
        odds_ratio = unname(ft$estimate)
    )
}))

resultado[, padj := p.adjust(pvalue, method = "BH")]


annot <- fread("/home/ttamaki/New/results/interproscan/genbank_interproscan_annotations.tsv", header = TRUE)

cluster_annot <- merge(
  member_id,
  annot[, .(
    Acession,
    `InterPro Acession`,
    `InterPro Description`
  )],
  by.x = "gene",
  by.y = "Acession"
)

r2 <- merge(
    resultado,
    unique(
        cluster_annot[ , .(
            cluster_n,
            `InterPro Acession`,
            `InterPro Description`)
            ]),
    by.x = "cluster_id",
    by.y = "cluster_n"
)
r2 <- r2[,.(
    cluster_id,
    `InterPro Acession`,
    `InterPro Description`,
    `odds_ratio`,
    `padj`
)]
fwrite(r2, "/home/ttamaki/New/results/teste_fisher_r7.6_c0.16_t4cp.tsv", sep = "\t")
r3 <- r2[padj<0.05][order(odds_ratio, decreasing = TRUE)]
fwrite(r3, "/home/ttamaki/New/results/teste_fisher_r7.6_c0.16_t4cp_a0.05.tsv", sep = "\t")
#maiores20 <- resultado[padj<0.05][odds_ratio!=Inf][order(odds_ratio, decreasing = TRUE)][1:20]
#maiores20
#olhar <- cluster_annot[cluster_n %in% maiores20$cluster_id]
#View(unique(olhar$`InterPro Description`))
#menores20 <- resultado[order(padj)][1:20]
#olhar <- cluster_annot[cluster_n %in% menores20$cluster_id]    
#olhar


