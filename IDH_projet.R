library(STRINGdb)
library(tidyverse)
packageVersion("STRINGdb")


sdb = STRINGdb$new(version='10', species = 6183, score_threshold=400, input_directory='repo_data')
g = sdb$get_graph()
sp = sdb$get_proteins()
annot = sdb$get_annotations()



library(neo4r)
packageVersion("neo4r")

con = neo4j_api$new(
  url="http://localhost:7474",
  user="neo4j",
  password ="compolfecs"  #Pas terrible ça, c'est un de mes mdp mais bon pas trop le choix
)

#Querys for proteins and GOterms
query_prots = paste0("MATCH (n:Protein) RETURN n")
prots = query_prots %>% call_neo4j(con)

query_terms = paste0("MATCH (m:GOTerm) RETURN m")
termsGO = query_terms %>% call_neo4j(con)

# Example
# Get proteins from a specific GO term
goterm = 'GO:0008150'
query = paste0( "MATCH (g:GOTerm {acc: '",goterm,"'})-[*] -> (p:Protein) return p")
proteins = query %>% call_neo4j(con)


# On veut pour chaque GOTerm, la liste des gènes annotés sur ce GO
# De façon directe ou implicite (direct + les prots des GO termes descendants)
query_annot = paste0("MATCH (n:Protein)<-[r:GOAnnotates]-(m:GOTerm) RETURN m, r, n")
prot_annnot_term = query_annot %>% call_neo4j(con)
#Les GO termes annotés sur des protéines
GO_terms_annotations = prot_annnot_term$m$acc
#Les GO termes uniques annotés sur des protéines
GO_terms_annotations_uniq = unique(prot_annnot_term$m$acc)

#On écrit le fichier des protéines directement et implicitement annotées par les GOtermes
fileConn<-file("Prot_names_imp.txt")
writeLines(c("# format: sets", "# comment : Biological sets of proteins annotated with GO terms implictly (with more specific descending GO terms included)"), fileConn)
close(fileConn)
cat("#Species : ", file = "Prot_names_imp.txt","\t", append =TRUE)
cat(sdb$species, file = "Prot_names_imp.txt", "\n", append =TRUE)
for (go_term in GO_terms_annotations_uniq) {
  goterm=go_term
  query_dir_imp = paste0( "MATCH (g:GOTerm {acc: '",goterm,"'})-[*] -> (p:Protein) return g.term_type, g.name, p.protein_external_id")
  proteins_dir_imp = query_dir_imp %>% call_neo4j(con)
  protein_names_dir_imp = proteins_dir_imp$p$value
  protein_names_dir_imp = unique(protein_names_dir_imp)
  go_term_type = proteins_dir_imp$g.term_type
  go_type = go_term_type$value[1]
  go_name_type = proteins_dir_imp$g.name
  go_name = go_name_type$value[1]
  if (length(protein_names_dir_imp)!=0) {
    cat(paste0(goterm, "\t", go_type, go_name), file ="Prot_names_imp.txt", append = TRUE)
    cat(substr(protein_names_dir_imp,6,100), file = "Prot_names_imp.txt", append = TRUE)
    cat("\n", file = "Prot_names_imp.txt", append = TRUE)}
}

#On écrit le fichier des protéines directement annotées par les GOtermes
fileConn<-file("Prot_names_dir.txt")
writeLines(c("# format: sets","# comment : Biological sets of proteins directly annotated with GOTerms"), fileConn)
close(fileConn)
cat("#Species : ", file = "Prot_names_dir.txt","\t", append =TRUE)
cat(sdb$species, file = "Prot_names_dir.txt", "\n", append =TRUE)
for (go_term in GO_terms_annotations_uniq) {
  goterm=go_term
  query_dir = paste0( "MATCH (g:GOTerm {acc: '",goterm,"'})--> (p:Protein) return g.term_type, g.name, p.protein_external_id") 
  proteins_dir = query_dir %>% call_neo4j(con)
  protein_names_dir = proteins_dir$p$value
  protein_names_dir = unique(protein_names_dir)
  go_term_type = proteins_dir$g.term_type
  go_type = go_term_type$value[1]
  go_name_type = proteins_dir$g.name
  go_name = go_name_type$value[1]
  if (length(protein_names_dir)!=0) {
    cat(paste0(goterm, "\t", go_type, go_name), file ="Prot_names_dir.txt", append = TRUE)
    cat(substr(protein_names_dir,6,100), file = "Prot_names_dir.txt", append = TRUE)
    cat("\n", file = "Prot_names_dir.txt", append = TRUE)
}}

  



# On écrit le fichier Résumé de l'annotation
fileConn=file("Resume_GO.txt")
writeLines(c("# format: Resume"), fileConn)
close(fileConn)
cat("#Species : ", file = "Resume_GO.txt", "\t", append =TRUE)
cat(sdb$species, file = "Resume_GO.txt", "\n", append =TRUE)
cat("Nombre d'annotations totales", file = "Resume_GO.txt", "\t", append =TRUE)
nb_annot = nrow(annot)
cat(nb_annot, file = "Resume_GO.txt", "\n", append =TRUE)
cat("Nombre d'annotations GO", file = "Resume_GO.txt", "\t", append =TRUE)
GO_annots = length(which(startsWith(annot$term_id,'GO'))==TRUE)
cat(GO_annots, file = "Resume_GO.txt", "\n", append =TRUE)
cat("Nombre de protéines uniques annotées", file = "Resume_GO.txt", "\t", append =TRUE)
uniq_annot_prot = length(unique(prot_annnot_term$n$name))
cat(uniq_annot_prot, file = "Resume_GO.txt", "\n", append =TRUE)
cat("Nombre de termes GO annotés sur une protéine", file = "Resume_GO.txt", "\t", append =TRUE)
nb_prots_GO = length(prot_annnot_term$m$acc)
cat(nb_prots_GO, file = "Resume_GO.txt", "\n", append =TRUE)
cat("Nombre de termes par gène", file = "Resume_GO.txt", "\t", append =TRUE)
terms_per_gene = nb_prots_GO/uniq_annot_prot
cat(terms_per_gene, file = "Resume_GO.txt", "\n", append =TRUE)
cat("Nombre de gènes annotés par terme", file = "Resume_GO.txt", "\t", append =TRUE)
genes_per_term = uniq_annot_prot /nb_prots_GO
cat(genes_per_term, file = "Resume_GO.txt", "\n", append = TRUE)
cat("Nombre de protéines totales", file = "Resume_GO.txt", "\t", append =TRUE)
prots_all = length(unique(prots$n$name))
cat(prots_all, file = "Resume_GO.txt", "\n", append =TRUE)
cat("Nombre de termes GO uniques annotés sur une protéine", file = "Resume_GO.txt", "\t", append =TRUE)
uniq = length(unique(prot_annnot_term$m$acc))
cat(uniq, file = "Resume_GO.txt", "\n", append = TRUE)



