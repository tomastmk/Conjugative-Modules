library(data.table)
library(tidyr)


# Módulos do CONJScan
CONJScan <- fread("/home/ttamaki/Conj/data/CONJScan_modules.tsv", header = FALSE, sep = "\t")
head(CONJScan)
CONJScan_exploded <- CONJScan %>% 
  separate_longer_delim(V2, delim = ",")


CONJScan_genes <- paste0(CONJScan_exploded$V1,"_",CONJScan_exploded$V2)

# Comunidades
member_id <- fread("/home/ttamaki/New/results/member_id.tsv", header = TRUE)

head(member_id)
max(member_id$cluster_n)

clusters_communities_path <- "/home/ttamaki/New/results/var_com/jaccard_com_r5.0_c0.016.tsv" 
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
communities[, conjugative := any(gene %in% CONJScan_genes), by = community_id]
conjugative_communities <- communities[conjugative == TRUE]
conjugative_communities[, conjugative := NULL]

### Adiciona o id do módulo
conjugative_communities[
    ,
    module_id := .GRP,
    by = .(community_id, plasmid)
]

### Tira módulos que possuem apenas um gene
conjugative_communities <- conjugative_communities[
    conjugative_communities[, .I[.N > 1], by = module_id]$V1
]

conjugative_communities[, in_conjscan := gene %in% CONJScan_genes]

module_info <- conjugative_communities[
    ,
    .(conj = any(in_conjscan)),
    by = .(community_id, module_id)
]


module_cluster <- unique(
    conjugative_communities[
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
                conj,
                presente = cluster_id == cl
            )
        ]
    )

    presente <- presente[
        ,
        .(presente = any(presente)),
        by = .(community_id, module_id, conj)
    ]

    ft <- fisher.test(table(presente$conj, presente$presente),
                      alternative = "greater")

    data.table(
        cluster_id = cl,
        pvalue = ft$p.value,
        odds_ratio = unname(ft$estimate)
    )
}))
resultado[, padj := p.adjust(pvalue, method = "BH")]
fwrite(resultado, "/home/ttamaki/New/results/teste_fisher.tsv", sep = "\t")



resultado <- fread("/home/ttamaki/New/results/teste_fisher.tsv", header = TRUE)
annot <- fread("/home/ttamaki/New/results/interproscan/genbank_interproscan_annotations.tsv", header = TRUE)
member_id <- fread("/home/ttamaki/New/results/member_id.tsv", header = TRUE)
unique(member_id[cluster_n==125254]$gene %in% CONJScan_genes)




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
fwrite(r2, "/home/ttamaki/New/results/teste_fisher_annot.tsv", sep = "\t")
maiores20 <- resultado[padj<0.05][odds_ratio!=Inf][order(odds_ratio, decreasing = TRUE)][1:20]
maiores20
olhar <- cluster_annot[cluster_n %in% maiores20$cluster_id]
View(unique(olhar$`InterPro Description`))
menores20 <- resultado[order(padj)][1:20]
olhar <- cluster_annot[cluster_n %in% menores20$cluster_id]
olhar
