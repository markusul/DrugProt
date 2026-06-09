## app.R ##
library(shinydashboard)
library(plotly)
library(networkD3)
library(grid)
library(readxl)
library(shinyWidgets)
library(pryr)

#source("dbAccess.R")
#db   <- dp_connect()
#meta <- dp_metadata(db)

source("dbAccess_parquet.R")
db <- dp_connect("parquet", in_memory_network = FALSE, memory_limit = "3GB")
meta <- dp_metadata(db)

replace_drug_ids <- function(x) {
  # 1. Split the string by backticks and underscores
  # E.g., "`drug_#123`" becomes c("", "drug", "#123", "")
  # E.g., "`drug_EV`"   becomes c("", "drug", "EV", "")
  ids <- strsplit(x, "_|`")[[1]]
  
  # 2. Filter out empty strings and the "drug" prefix
  ids <- ids[ids != "" & ids != "drug" & ids != ":"]
  
  # If there's nothing left, return the original input
  if(length(ids) == 0) return(x)
  
  # 3. Process each remaining part
  names <- sapply(ids, function(id) {
    # Check if this part is a "number" (ID). 
    # We remove the '#' for the numeric check to handle both '123' and '#123'.
    is_numeric_id <- grepl("#", id) || !is.na(as.numeric(gsub("#", "", id)))
    
    if (is_numeric_id) {
      # It's an ID: Try to get the name from drug_lookup
      val <- drug_lookup[id]
      # If it's found in the table, return the name; otherwise return the ID itself
      if (!is.na(val)) return(unname(val))
      return(id)
    } else {
      # It's already a name (like "EV", "LA", "MK"): Return it as is
      return(id)
    }
  })
  
  # 4. Join parts back with ":" (useful for drug combinations)
  paste(names, collapse = ":")
}

prot_names_short <- meta$prot_names_short
treatment        <- meta$treatment
treatNames       <- meta$treatNames
expTimes         <- meta$expTimes
drugOrder        <- meta$drugOrder
nProtein         <- meta$nProtein
nTreatment       <- meta$nTreatment
nDrugs           <- length(drugOrder)

tmp <- dbGetQuery(db, "SELECT id, name FROM drug_lookup")
drug_lookup <- setNames(tmp$name, tmp$id)

drugOrder <- sapply(drugOrder, replace_drug_ids)
t_choice  <- paste0(expTimes, "h")


ui <- dashboardPage(
  title = "Drug-Prot",
  
  dashboardHeader(title = tags$a(href='https://ulme.shinyapps.io/DrugProt/',
                                 tags$img(src='drugprot-logo.svg', height = '50'))
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("About", tabName = "About", icon = icon("circle-info")),
      menuItem("Settings", tabName = "Settings", icon = icon("cog")),
      menuItem("Drug Effects", tabName = "DrugEffects", icon = icon("dashboard")),
      menuItem("Protein Network", tabName = "ProteinNetwork", icon = icon("th"))
    )
  ),
  dashboardBody(
    withMathJax(),
    tabItems(
      tabItem(tabName = "About",
              h2("Drug-Prot: statistical inference of drug effects and protein dependencies"),
              fluidRow(
                box(title = "What Drug-Prot does", status = "primary", solidHeader = TRUE, width = 12,
                    p("Drug-Prot is an interactive query system built on a large-scale perturbation proteomics dataset of 18 breast cancer cell lines, ",
                      "treated with 63 single drugs and 59 drug combinations, with protein expression measured at 6, 24, and 48 hours after treatment ",
                      "(Sun et al., 2025). For a user-defined set of proteins, Drug-Prot reports pre-computed statistical evidence for two kinds of relationship:"),
                    tags$ul(
                      tags$li(strong("Drug effects on proteins:"), " corrected p-values for the effect of each single drug and each drug pair on the selected proteins, at each time point. ",
                              "These are reported on the ", strong("Drug Effects"), " tab."),
                      tags$li(strong("Temporal dependencies between proteins:"), " a directed network in which an edge from one protein to another indicates that the earlier protein's differential expression ",
                              "is significantly associated with the later protein's, after adjusting for residual drug effects. These are reported on the ", strong("Protein Network"), " tab.")
                    ),
                    p("All evidence is pre-computed (approximately 62 million p-values across all 5,392 measured proteins, 122 treatments, and three time points), ",
                      "so queries return instantly and you never need to download the underlying dataset."),
                    p("Throughout, proteins are referred to by their ", strong("HGNC gene symbol"), " (e.g. RAB7A, LMNA, MAP2K1). ",
                      "Where several proteins were quantified together as one group, their symbols are joined with a slash (e.g. CALM1/CALM2/CALM3).")
                )
              ),
              fluidRow(
                box(title = "How to interpret the results", status = "warning", solidHeader = TRUE, width = 12,
                    p(strong("The two kinds of relationship are not interpreted the same way.")),
                    tags$ul(
                      tags$li(strong("Drug \u2192 protein effects can be read causally."),
                              " Because drug administration is externally controlled and unconfounded, a significant effect can be interpreted as the interventional change in protein expression that the drug (or drug pair) would induce."),
                      tags$li(strong("Protein \u2192 protein edges are directional, but not necessarily causal."),
                              " Edge direction reflects the temporal ordering of the measurements (earlier \u2192 later), not a verified mechanism. ",
                              "Hidden, unmeasured biological processes acting between time points may influence both endpoints, so these edges are best treated as ",
                              em("hypothesis-generating"), " rather than mechanistic.")
                    ),
                    p("The drug effects (at 6 hours) and the protein dependencies (at 24 and 48 hours) are estimated from two high-dimensional linear models, sketched below for a single protein:"),
                    p("$$ y^{6} = \\sum_{j} (\\alpha^{6_{0}}_{j} + \\alpha^{6}_{j} D_{j}) + \\sum_{j \\neq k} (\\beta^{6_{0}}_{jk} + \\beta^{6}_{jk} D_{j} D_{k}) + \\varepsilon^{6} $$"),
                    p("$$ y^{t} = \\sum_{j} (\\alpha^{t_{0}}_{j} + \\alpha^{t}_{j} D_{j}) + \\sum_{j \\neq k} (\\beta^{t_{0}}_{jk} + \\beta^{t}_{jk} D_{j} D_{k}) + (Y^{t-})^{\\top} \\gamma^{t} + \\varepsilon^{t}, \\quad t \\in \\{24, 48\\} $$"),
                    p("where ", tags$code("y\u1d57"), " is the differential expression of the protein at time ", tags$code("t"), " (relative to its untreated baseline); ",
                      tags$code("D\u2c7c"), " is the administered concentration of drug ", tags$code("j"), "; ",
                      tags$code("Y\u1d57\u207b"), " is the vector of differential expressions of all measured proteins at the preceding time point; and ", tags$code("\u03b5"), " is noise. ",
                      "The ", tags$code("\u03b1\u2070"), ", ", tags$code("\u03b2\u2070"), " terms are treatment intercepts; ", tags$code("\u03b1"), ", ", tags$code("\u03b2"),
                      " capture single-drug and drug-interaction effects; and ", tags$code("\u03b3"), " captures the temporal protein-to-protein dependencies."),
                    p("Parameters are estimated with the de-sparsified Lasso, with group p-values for the drug terms and p-values for the protein terms. ",
                      "This is a simplified sketch: the full models, the precise definition of the aggregated differential expression ", tags$code("Y"),
                      " (which handles the unpaired measurements across time points), and all assumptions are given in the accompanying paper.")
                )
              ),
              fluidRow(
                box(title = "How a query works", status = "primary", solidHeader = TRUE, width = 6,
                    tags$ol(
                      tags$li("On the ", strong("Settings"), " tab, choose a set of proteins of interest \u2014 type/select them, click an ",
                              em("Example"), " button to load either the two most IC50-predictive proteins (LMNA and RAB7A) or the full 26-protein set from the paper, or upload a .txt file with one protein name per line."),
                      tags$li("Choose a significance level and, if you wish, change the multiple-testing correction methods (applied separately to drug and protein effects)."),
                      tags$li("Read off drug effects on the ", strong("Drug Effects"), " tab and the dependency network on the ", strong("Protein Network"), " tab."),
                      tags$li("Download any of the p-value tables (CSV) or the networks (interactive HTML) from the respective tabs.")
                    ),
                    p(em("Note:"), " for a queried set, both the parents and children of each protein are searched across the whole proteome, so the returned network can extend well beyond the proteins you selected.")
                ),
                box(title = "Access & links", status = "primary", solidHeader = TRUE, width = 6,
                    p(strong("Paper: "), a(href = "#", "[link to be added]", target = "_blank")),
                    p(strong("Web application: "), a(href = "https://ulme.shinyapps.io/DrugProt/", "ulme.shinyapps.io/DrugProt", target = "_blank")),
                    p(strong("Software source code: "), a(href = "https://github.com/markusul/DrugProt", "github.com/markusul/DrugProt", target = "_blank")),
                    p(strong("Paper code (p-value computation): "), a(href = "https://github.com/markusul/SDForest-Paper", "github.com/markusul/SDForest-Paper", target = "_blank")),
                    p(strong("Underlying dataset: "), "Sun et al. (2025), ",
                      a(href = "https://doi.org/10.1101/2025.02.07.637070", "doi.org/10.1101/2025.02.07.637070", target = "_blank"))
                )
              ),
              fluidRow(
                box(title = "Cite Drug-Prot", status = "primary", solidHeader = TRUE, width = 12,
                    p("If you use Drug-Prot, please cite:"),
                    p(em("Ulmer, M., Sun, R., Qian, L., Aebersold, R., Guo, T., and B\u00fchlmann, P. (2026). ",
                         "Drug-Prot: A query system for statistical inference of drug effects and interactions in dynamic proteomic networks.")),
                    p("Please also cite the underlying dataset:"),
                    p(em("Sun, R., Qian, L., Li, Y., et al. (2025). A perturbation proteomics-based foundation model for virtual cell construction. ",
                         "bioRxiv. https://doi.org/10.1101/2025.02.07.637070"))
                )
              ),
              fluidRow(
                box(title = "Methods & references", status = "primary", solidHeader = TRUE, width = 12, collapsible = TRUE, collapsed = TRUE,
                    p("The statistical methods underlying Drug-Prot are described in detail in the accompanying paper. Key references:"),
                    tags$ul(
                      tags$li("Zhang, C.-H. and Zhang, S. S. (2014). Confidence intervals for low dimensional parameters in high dimensional linear models. ",
                              em("J. R. Stat. Soc. B"), " 76(1):217\u2013242. (de-sparsified Lasso)"),
                      tags$li("van de Geer, S., B\u00fchlmann, P., Ritov, Y., and Dezeure, R. (2014). On asymptotically optimal confidence regions and tests for high-dimensional models. ",
                              em("Ann. Statist."), " 42(3):1166\u20131202."),
                      tags$li("Dezeure, R., B\u00fchlmann, P., Meier, L., and Meinshausen, N. (2015). High-dimensional inference: confidence intervals, p-values and R-software hdi. ",
                              em("Statist. Sci."), " 30(4):533\u2013558."),
                      tags$li("B\u00fchlmann, P. (2013). Statistical significance in high-dimensional linear models. ",
                              em("Bernoulli"), " 19(4):1212\u20131242. (group p-values)"),
                      tags$li("Holm, S. (1979). A simple sequentially rejective multiple test procedure. ",
                              em("Scand. J. Statist."), " 6(2):65\u201370. (FWER control for drug effects)"),
                      tags$li("Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate. ",
                              em("J. R. Stat. Soc. B"), " 57(1):289\u2013300. (FDR control for protein effects)")
                    )
                )
              )
      ),
      tabItem(tabName = "Settings",
              h2("Settings"),
              fluidRow(
                box(title = "Select Proteins of Interest", status = "primary", solidHeader = TRUE,
                    pickerInput("protSet", "Select Protein Set", 
                                choices = unname(prot_names_short), 
                                multiple = TRUE, 
                                options = list("live-search"=TRUE)),
                    helpText("Start typing to search the 5,392 measured proteins, and select one or more to query. Proteins are named by their HGNC gene symbol (e.g. RAB7A, LMNA, MAP2K1); proteins quantified together as a group are joined with a slash (e.g. CALM1/CALM2/CALM3). The analysis is restricted to the selected set, which reduces the multiple-testing burden and increases power."),
                    actionButton("preSelectedTwo", "Example: LMNA & RAB7A"),
                    actionButton("preSelected", "Example: 26 IC50 proteins"),
                    actionButton("clear", "Clear Selection"), 
                    helpText("\"Example: LMNA & RAB7A\" loads the two most IC50-predictive proteins from the paper \u2014 a small, interpretable network. \"Example: 26 IC50 proteins\" loads the full predictive set from the paper. \"Clear Selection\" empties the current set."),
                    fileInput("file", "Upload .txt File with Protein Names (one per line)", accept = c(".txt")),
                    helpText("Alternatively, upload a plain-text file with one protein name per line (e.g. a curated pathway or complex). Names must match those in the dataset."),
                    h3("p-value Adjustments"),
                    numericInput("alpha", "Significance Level", value = 0.05, min = 0, max = 1, step = 0.0001),
                    helpText("The threshold below which a corrected p-value is called significant (default 0.05)."),
                    selectInput("corectionDrug", "Correction Method for Drug Effects", choices = p.adjust.methods, selected = 'holm'),
                    helpText("Multiple-testing correction applied to drug effects across the selected proteins and time points. Holm controls the family-wise error rate (default)."),
                    selectInput("corectionProtein", "Correction Method for Protein Effects", choices = p.adjust.methods, selected = 'BH'),
                    helpText("Correction applied separately to the protein-to-protein dependencies. Benjamini-Hochberg controls the false discovery rate (default), which is better suited to the much larger number of protein-protein tests."),
                    width = 6), 
                box(title = "Download P-values for selected Proteins", status = "primary", solidHeader = TRUE,
                  downloadButton("downloadPvalsDrug", "Download P-values of Drug Effects"), 
                  downloadButton("downloadPvalsProtein", "Download P-values of Protein Effects"), width = 6),
                box(title = "Selected Proteins", status = "primary", solidHeader = TRUE,
                    tableOutput("selectedTable"), width = 3),
                box(title = "Summary", status = "primary", solidHeader = TRUE,
                    tableOutput("summaryTable"), width = 3)
              )
      ),
      tabItem(tabName = "DrugEffects", 
              h2("Drug Effects"),
              checkboxGroupButtons("t", "Select time after drug administration to analyze", t_choice, t_choice, checkIcon = list(
                yes = icon("square-check"),
                no = icon("square"))),
              helpText("Select one or more post-treatment time points (6, 24, 48 hours). Reported p-values are the minimum over the selected time points and the queried proteins."),
              fluidRow(
                box(title = "How to read the heatmap", status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
                    tags$ul(
                      tags$li(strong("Diagonal cells"), " show the evidence that a ", strong("single drug"), " affects the selected protein set."),
                      tags$li(strong("Off-diagonal cells"), " show the evidence for an ", strong("interaction"), " between a pair of drugs on the selected protein set."),
                      tags$li(strong("Grey scale"), " encodes the corrected p-value (lighter = smaller p-value = stronger evidence)."),
                      tags$li(strong("Black cells"), " mark drug pairs for which no experimental data exist (shown with a sentinel value of 2); they are not tested.")
                    )
                )
              ),
              fluidRow(
                box(title = "Drug Effect Heatmap", status = "primary", solidHeader = TRUE,
                    plotlyOutput("plotDrugEffects", height = 800), width = 6,
                    selectInput("selTreat", "Effects of treatment", choices = treatNames),
                    helpText("We show the average treatment effect for the concentrations in the data at the time point of highest significance."),
                    helpText("The format is: effect (time point p-value)."),
                    tableOutput("DrugEffects")
                    ),
                box(title = "Significant Single Drugs", status = "primary", solidHeader = TRUE,
                    tableOutput("singleDrugs"), width = 3),
                box(title = "Significant Drug Interactions", status = "primary", solidHeader = TRUE,
                    tableOutput("interactions"), width = 3)
              )
      ),
      tabItem(tabName = "ProteinNetwork", 
              h2("Protein Network"),
              fluidRow(
                box(title = "How to read the network", status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
                    tags$ul(
                      tags$li("A directed edge from one protein to another indicates a ", strong("temporal dependency"),
                              ": the earlier protein's differential expression is significantly associated with the later protein's, after adjusting for residual drug effects."),
                      tags$li("In the ", strong("Summary Graph"), ", nodes are coloured by group: the proteins you queried (", em("Selected"),
                              ") and the further proteins drawn in because they are significantly linked to your set (", em("Connected"), ")."),
                      tags$li("In the ", strong("Temporal Graph"), ", each node is a (protein, time-point) pair, coloured by time point (6h, 24h, 48h), so the same protein can appear at several times."),
                      tags$li("In the ", strong("Temporal Graph"), ", each connection is coloured according to the sign of the estimated effect", tags$code("y\u1d57"), ". Negative effects are ", tags$span(style = "color: red; font-weight: bold;", "red"), ", while positive effects are ", tags$span(style = "color: blue; font-weight: bold;", "blue"), "."),
                      tags$li("The network is searched across the ", strong("whole proteome"), ", so it routinely extends well beyond the proteins you selected."),
                      tags$li(em("Reminder:"), " edge direction reflects temporal ordering, not a verified causal mechanism, and may be subject to unmeasured confounding.")
                    )
                )
              ),
              fluidRow(
                box(title = "Summary Graph", status = "primary", solidHeader = TRUE,
                    forceNetworkOutput("SummaryGraph", height = 800), width = 12)
              ),
              fluidRow(
                box(title = "Temporal Graph", status = "primary", solidHeader = TRUE,
                    forceNetworkOutput("TemporalGraph", height = 800), width = 12)
              ),
              fluidRow(
                box(title = "Download Protein Interaction Network", status = "primary", solidHeader = TRUE,
                    downloadButton("downloadSummary", "Download Summary Graph"), 
                    downloadButton("downloadTemporal", "Download Temporal Graph"), width = 4),
                box(title = "Relevant Proteins", status = "primary", solidHeader = TRUE,
                    fluidRow(column(tableOutput("numProteinEffects"), width = 6), 
                             column(tableOutput("ProteinEffects"), width = 6)), width = 8)
              )
      )
    )

  )
)

server <- function(input, output) {
  # control panel
  print("Start Server")

  observeEvent(input$preSelectedTwo, {
    updatePickerInput(session = getDefaultReactiveDomain(), inputId = "protSet", 
                      selected = c("LMNA", "RAB7A"))
  })
  observeEvent(input$preSelected, {
    most_imp_26 <- c("VDAC1", "CALM1/CALM2/CALM3", "IQGAP1", "RPS25", "PCBP2", "DDX39B",
                     "GMPS", "MAP2K1", "KIF5B", "RPLP0", "CYCS", "SNRNP70", "PTMA", "LMNA",
                     "SFPQ", "ACLY", "RALY", "RPSA", "RPL7", "SF3B3", "RAB7A", "SUPT16H",
                     "MYL6", "RAP1B", "AKAP13", "HSP90B1")
    updatePickerInput(session = getDefaultReactiveDomain(), inputId = "protSet", 
                      selected = prot_names_short[prot_names_short %in% most_imp_26])
  })
  observeEvent(input$clear, {
    updatePickerInput(session = getDefaultReactiveDomain(), inputId = "protSet", selected = character(0))
  })
  observeEvent(input$file, {
    req(input$file)
    prot_upload <- readLines(input$file$datapath)
    prot_upload <- prot_upload[prot_upload %in% prot_names_short]
    updatePickerInput(session = getDefaultReactiveDomain(), inputId = "protSet", selected = prot_upload)
  })
  
  protSet <- reactive({
    if (is.null(input$protSet))
      c()
    else
      input$protSet
  })
  protSet_d <- protSet %>% debounce(1300)
  
  P_selection <- reactive({
    print(protSet_d())
    sel <- which(prot_names_short %in% protSet_d())
    if(length(sel) == 0) return(NULL)
    sel
  })

  
  adjusted_pvecs <- reactive({
    if(is.null(P_selection())) return(NULL)
    print("Fetching and adjusting p-values from DB...")
    
    selPvecs <- dp_drug_selPvecs(db, P_selection(), nTreatment, length(expTimes))
    array(p.adjust(selPvecs, method = input$corectionDrug), dim = dim(selPvecs))
  })
  
  pvec <- reactive({
    if(is.null(adjusted_pvecs()) || is.null(input$t)) return(NULL)
    t_selection <- t_choice %in% input$t

    selPvecs <- adjusted_pvecs()
    
    pvec <- apply(selPvecs, 1, function(p) min(p[t_selection, ]))
    names(pvec) <- sapply(treatment, replace_drug_ids)
    pvec
  })
  
  min_loc <- reactive({
    if(is.null(adjusted_pvecs()) || is.null(input$t)) return(NULL)
    
    t_selection <- t_choice %in% input$t
    selPvecs <- adjusted_pvecs()
    
    min_sel <- apply(selPvecs, c(1, 3), function(p) which.min(p[t_selection]))

    apply(selPvecs, c(1, 3), function(p) which.min(p[t_selection]))
  })
  
  Links_all <- reactive({
    print("links")
    if(is.null(P_selection())) return(NULL)
    selection <- P_selection()
    
    # select relevant p-values
    print(mem_used())
    Pval_sel <- dp_net_pvals(db, selection, input$alpha, nProtein, length(expTimes) - 1)
    
    print("correction")
    # p-value correction
    Tpval <- lapply(Pval_sel, function(links) links$pvalue)
    nPpertime <- sapply(Tpval, length)
    
    # number of p values to correct for
    nPvalues <- (length(expTimes) - 1)*length(selection)*(2*nProtein - length(selection))
    
    Tpval <- p.adjust(unlist(Tpval), method = input$corectionProtein, 
                                  n = nPvalues)
    
    # Dynamically assign the adjusted p-values back to their respective time points
    start_idx <- 1
    for (i in seq_along(Pval_sel)) {
      end_idx <- start_idx + nPpertime[i] - 1
      
      Pval_sel[[i]][, "pvalue"] <- Tpval[start_idx:end_idx]
      
      # Update the starting index for the next time point
      start_idx <- end_idx + 1
    }
    
    # convert p-values to links
    Links_all <- lapply(Pval_sel, function(links) links[links$pvalue < input$alpha, ])
    
    if(sum(sapply(Links_all, nrow)) == 0) return(NULL)
    print(mem_used())
    Links_all
  })

  SummGraph <- reactive({
    print("SumGraph")
    if(is.null(Links_all())) return(NULL)
    Links_sum <- do.call(rbind, Links_all())
    Links_sum[, "source"] <- Links_sum[, "source"] - 1
    Links_sum[, "target"] <- Links_sum[, "target"] - 1
    Links_sum[, "value"] <- 1
  
    rel.Nodes <- sort(unique(c(P_selection()-1, unlist(Links_sum[, c('source', 'target')]))))
    Nodes_sum <- data.frame(name = prot_names_short[rel.Nodes+1], group = "Connected", size = 1)
    Nodes_sum$group[rel.Nodes %in% (P_selection()-1)] <- "Selected"

    #reorganize link index
    for(i in 1:length(rel.Nodes)){
      Links_sum[, c('source', 'target')][Links_sum[, c('source', 'target')] == rel.Nodes[i]] <- i - 1
    } 
    list(Links_sum = Links_sum, Nodes_sum = Nodes_sum)
  })
  
  # temporal graph
  TempGraph <- reactive({
    print("TempGraph")
    # if(is.null(Links_all())) return(NULL)
    
    nT <- length(expTimes)
    
    # 1. Dynamically build the 'rel' list (nodes present at each time step)
    rel <- list()
    for(t in 1:nT) {
      if (t == 1) {
        # First time point: only sources
        rel[[t]] <- sort(unique(c(P_selection(), Links_all()[[t]][, "source"])))
        
      } else if (t == nT) {
        # Last time point: only targets
        rel[[t]] <- sort(unique(c(P_selection(), Links_all()[[t-1]][, "target"])))
      } else {
        # Middle time points: targets from the previous step, sources for the next step
        rel[[t]] <- sort(unique(c(P_selection(), Links_all()[[t-1]][, "target"], Links_all()[[t]][, "source"])))
      }
    }
    
    # 2. Calculate continuous index offsets for D3
    lenRel <- c(0, sapply(rel, length))
    
    # 3. Create node names and groups dynamically
    nodenames <- unlist(lapply(1:nT, function(t) paste(prot_names_short[rel[[t]]], expTimes[t], sep = '_')))
    nodegroups <- rep(paste0(expTimes, "h"), times = sapply(rel, length))
    
    # 4. Build the temporal links
    if(is.null(Links_all())){
      Links_temp <- data.frame(source = 0, target = 0, value = 1)
    } else {
      # Loop over the transitions (number of timepoints - 1)
      Links_temp <- lapply(1:(nT - 1), function(t) {
        links <- Links_all()[[t]]
        
        # Re-index sources for this transition to match the flattened node list
        for(i in seq_along(rel[[t]])){
          links[links[, 1] == rel[[t]][i], 1] <- i - 1 + sum(lenRel[1:t])
        }
        # Re-index targets for this transition
        for(i in seq_along(rel[[t+1]])){
          links[links[, 2] == rel[[t+1]][i], 2] <- i - 1 + sum(lenRel[1:(t+1)])
        }
        return(links)
      })
      Links_temp <- do.call(rbind, Links_temp)
      Links_temp$value <- 1
    }
    
    # 5. Assemble final Nodes dataframe
    Nodes_temp <- data.frame(name = nodenames, group = nodegroups, size = 0.3)
    Nodes_temp$radius <- as.numeric(unlist(rel))
    
    
    #Pcoef  <- unlist(lapply(1:(nT-1), function(t){
    #  target <- unique(Links_all()[[t]][, "target"])
    #  unlist(lapply(target, function(protein){
    #    bhat <- dp_protein_bhat(db, protein, time_idx = t + 1, nProtein)
    #    unname(bhat[Links_all()[[t]][Links_all()[[t]][, "target"] == protein, "source"]])
    #    
    #  }))
    #}))
    
    Pcoef <- unlist(lapply(1:(nT-1), function(t) {
      L <- Links_all()[[t]]
      if (nrow(L) == 0) return(numeric(0))
      ti <- t + 1
      targets <- unique(L[, "target"])
      # ONE query for all targets in this transition:
      coef_df <- dp_protein_coef_bulk(db, targets, ti)   # cols: target, source, coef
      # look up each edge's coef by (target, source); absent -> 0
      key  <- paste(L[, "target"], L[, "source"], sep = "_")
      cmap <- setNames(coef_df$coef, paste(coef_df$target, coef_df$source, sep = "_"))
      out  <- cmap[key]
      out[is.na(out)] <- 0
      unname(out)
    }))

    direction <- rep("darkgrey", length(Pcoef))
    direction[Pcoef == 1] <- "blue"
    direction[Pcoef == 0] <- "red"
    
    print(mem_used())
    list(Links_temp = Links_temp, Nodes_temp = Nodes_temp, direction = direction)
  })

  output$selectedTable <- renderTable({
    if(length(P_selection()) == 0) return(NULL)
    data.frame(Proteins = prot_names_short[P_selection()], stringsAsFactors = FALSE)
  })

  output$summaryTable <- renderUI({
    if(length(P_selection()) == 0) return(NULL)
    sigEff <- pvec() < input$alpha

    P_sel_inNodes <- which(SummGraph()$Nodes_sum$name %in% prot_names_short[P_selection()]) - 1
    nParents <- sum(!SummGraph()$Links_sum$source %in% P_sel_inNodes)
    nChildren <- sum(!SummGraph()$Links_sum$target %in% P_sel_inNodes)
    
    drug_info <- c(
              paste0("Number of significant single drugs: ",
                     sum(sigEff[1:nDrugs])),
              paste0("Number of significant interactions: ",
                     sum(sigEff[(nDrugs + 1):length(sigEff)]))
    )

    protein_info <- c(
      paste0("Number of connections in protein network: ",
             nrow(SummGraph()$Links_sum)),
      paste0("Number of parents outside selected proteins: ",
             nParents),
      paste0("Number of children outside selected proteins: ",
             nChildren)
    )

    rows_drug <- paste(
      sapply(drug_info, function(line) {
        paste0("<tr><td>", line, "</td></tr>")
      }),
      collapse = ""
    )

    rows_protein <- paste(
      sapply(protein_info, function(line) {
        paste0("<tr><td>", line, "</td></tr>")
      }),
      collapse = ""
    )

    htmlTable <- paste0(
      "<table class=\"table table-bordered\">",
      "<tbody>",
      "<tr><th>Drug Effects</th></tr>",
      rows_drug,
      "<tr><th>Protein Network</th></tr>",
      rows_protein,
      "</tbody></table>"
    )

    HTML(htmlTable)
  })
  
  output$plotDrugEffects <- renderPlotly({
    if(length(P_selection()) == 0) return(NULL)
    if(is.null(pvec())) return(NULL)
    
    pMat <- matrix(NA, nrow = nDrugs, ncol = nDrugs)
    rownames(pMat) <- colnames(pMat) <- names(pvec())[1:nDrugs]
    for(l in names(pvec())){
      drugs <- strsplit(l, ":")[[1]]
      if(length(drugs) == 1) drugs <- c(drugs, drugs)
      pMat[drugs[1], drugs[2]] <- pvec()[l]
      pMat[drugs[2], drugs[1]] <- pvec()[l]
    }

    pMat <- as.matrix(pMat)
    pMat[is.na(pMat)] <- 2
    pMat <- pMat[drugOrder, drugOrder]

    ht <- plot_ly(z = pMat, x = colnames(pMat), y = colnames(pMat), 
                  type = "heatmap", colors = "Greys") %>%
                  layout(title = prot_names_short[P_selection()])
    ht
  })

  output$singleDrugs <- renderUI({
    if(length(P_selection()) == 0) return(NULL)
    if(is.null(pvec())) return(NULL)
    
    singleInd <- sapply(names(pvec()), function(l) length(strsplit(l, ":")[[1]]) == 1)
    pvec_single <- pvec()[singleInd]
    pvec_single <- pvec_single[order(pvec_single)]
    if(length(pvec_single) == 0) return(NULL)
    df <- data.frame(Drug = names(pvec_single), PValue = pvec_single, stringsAsFactors = FALSE)
    htmlTable <- paste0(
      '<table class="table table-bordered"><thead><tr><th>Drug</th><th>PValue</th></tr></thead><tbody>',
      paste(
        sapply(1:nrow(df), function(i) {
          pv <- df$PValue[i]
            color <- ifelse(pv < input$alpha, ' style="background-color: #add8e6;"', '')
          paste0('<tr><td', color, '>', df$Drug[i], '</td><td', color, '>', format(pv, digits = 4), '</td></tr>')
        }),
        collapse = ""
      ),
      '</tbody></table>'
    )
    HTML(htmlTable)
  })

  output$DrugEffects <- renderUI({
    print("drugEffects")
    if(length(P_selection()) == 0) return(NULL)
    if(is.null(pvec())) return(NULL)
    
    t_selection <- t_choice %in% input$t
    
    selTreat <- strsplit(input$selTreat, ":")[[1]]
    selTreat <- which(treatNames %in% c(selTreat, input$selTreat))

    df <- data.frame(Protein = unname(prot_names_short[P_selection()]), 
                     stringsAsFactors = FALSE)
    for(i in selTreat){
      dEff  <- dp_drug_dEff(db, i, nProtein, length(expTimes))
      DCoef <- dEff[P_selection(), t_selection]
      
      DCoef <- matrix(DCoef, nrow = length(P_selection()))
      DCoef <- unlist(lapply(1:length(P_selection()), function(t) DCoef[t, min_loc()[i, t]]))
      
      df <- cbind(df, DCoef)
    }
    names(df) <- c("Protein", names(treatNames[selTreat]))
    
    
    locations <- matrix(min_loc()[selTreat, ], nrow = length(selTreat))
    times <- t(apply(locations, 1:2, function(ii)t_choice[t_selection][ii]))
    
    
    seladjP <- matrix(0, nrow = nrow(times), ncol = ncol(times))
    for(i in 1:nrow(times)){
      for(j in 1:ncol(times)){
        seladjP[i, j] <- adjusted_pvecs()[selTreat[j], which(t_selection)[min_loc()[selTreat[j], i]], i]
      }
    }
    
    # Build the table header dynamically
    headers_html <- paste0(
      "<thead><tr>",
      paste(sapply(names(df), function(colname) paste0("<th>", colname, "</th>")), collapse = ""),
      "</tr></thead>"
    )
    
    # Build the table body dynamically
    body_html <- paste0(
      "<tbody>",
      paste(
        sapply(1:nrow(df), function(i) {
          # First column is always the Protein name (no color)
          row_html <- paste0('<tr><td>', df[i, 1], '</td>')
          
          # Loop through the remaining columns
          if (ncol(df) > 1) {
            cells_html <- sapply(2:ncol(df), function(j) {
              val_num <- df[i, j] # The coefficient (dEff)
              
              # Extract corresponding time and p-value
              # We use j - 1 because times and seladjP don't have the Protein column
              time_val <- times[i, j - 1]
              pval_num <- seladjP[i, j - 1]
              
              # Determine color based on positive/negative ONLY if p-value is significant
              if (!is.na(val_num) && !is.na(pval_num) && pval_num < input$alpha) {
                if (val_num < 0) {
                  color <- ' style="background-color: #ffcccc;"' # Light Red
                } else if (val_num > 0) {
                  color <- ' style="background-color: #cce5ff;"' # Light Blue
                } else {
                  color <- ''
                }
              } else {
                color <- '' # No color if it's NA or not significant
              }
              
              # Format the numbers nicely and merge them
              if (!is.na(val_num)) {
                formatted_coef <- format(val_num, digits = 4)
                # Using signif() for p-values keeps them clean, especially if very small
                formatted_pval <- signif(pval_num, digits = 3) 
                
                # Combine into: "coef (time pval)"
                display_str <- paste0(formatted_coef, " (", time_val, " ", formatted_pval, ")")
              } else {
                display_str <- ""
              }
              
              paste0('<td', color, '>', display_str, '</td>')
            })
            # Append the dEff cells to the row
            row_html <- paste0(row_html, paste(cells_html, collapse = ""))
          }
          
          row_html <- paste0(row_html, '</tr>')
          return(row_html)
        }),
        collapse = ""
      ),
      "</tbody>"
    )
    
    # Combine header and body into the final table
    htmlTable <- paste0('<table class="table table-bordered">', headers_html, body_html, '</table>')
    
    HTML(htmlTable)
  })
  
  output$interactions <- renderUI({
    if(length(P_selection()) == 0) return(NULL)
    if(is.null(pvec())) return(NULL)
    
    interactionInd <- sapply(names(pvec()), function(l) length(strsplit(l, ":")[[1]]) > 1)
    pvec_interaction <- pvec()[interactionInd]
    pvec_interaction <- pvec_interaction[order(pvec_interaction)]
    if(length(pvec_interaction) == 0) return(NULL)
    df <- data.frame(DrugInteraction = names(pvec_interaction), PValue = pvec_interaction, stringsAsFactors = FALSE)
    htmlTable <- paste0(
      '<table class="table table-bordered"><thead><tr><th>Drug Interaction</th><th>PValue</th></tr></thead><tbody>',
      paste(
        sapply(1:nrow(df), function(i) {
          pv <- df$PValue[i]
            color <- ifelse(pv < input$alpha, ' style="background-color: #add8e6;"', '')
          paste0('<tr><td', color, '>', df$DrugInteraction[i], '</td><td', color, '>', format(pv, digits = 4), '</td></tr>')
        }),
        collapse = ""
      ),
      '</tbody></table>'
    )
    HTML(htmlTable)
  })

  output$numProteinEffects <- renderTable({
    if(length(P_selection()) == 0) return(NULL)
    if(is.null(SummGraph())){
      df <- data.frame(Protein = prot_names_short[P_selection()], 
                       num.Parents = integer(1), num.Children = integer(1))
      return(df)
    }
    
    P_sel_inNodes <- which(SummGraph()$Nodes_sum$name %in% prot_names_short[P_selection()])
    df <- lapply(P_sel_inNodes - 1, function(i) {
      Protein <- SummGraph()$Nodes_sum$name[i + 1]
      Children <- sum(SummGraph()$Links_sum[, "source"] == i)#-1
      Parents <- sum(SummGraph()$Links_sum[, "target"] == i)#-1
      df <- data.frame(Protein = Protein, num.Parents = Parents, num.Children = Children)
      return(df)
    })
    df <- do.call(rbind, df)
    df
  })

  output$ProteinEffects <- renderTable({
    if(length(P_selection()) == 0) return(NULL)
    P_sel_inNodes <- which(SummGraph()$Nodes_sum$name %in% prot_names_short[P_selection()])
    fam <- lapply(P_sel_inNodes - 1, function(i) {
      Children <- SummGraph()$Nodes_sum$name[SummGraph()$Links_sum[, "target"][SummGraph()$Links_sum[, "source"] == i & SummGraph()$Links_sum[, "target"] != i] + 1]
      Parents <- SummGraph()$Nodes_sum$name[SummGraph()$Links_sum[, "source"][SummGraph()$Links_sum[, "target"] == i & SummGraph()$Links_sum[, "source"] != i] + 1]
      return(list(Children = Children, Parents = Parents))
    })

    # vector of all children and parents
    allChildren <- unique(unlist(sapply(fam, "[[", 1)))
    allParents <- unique(unlist(sapply(fam, "[[", 2)))

    # add dummies to smaller vectors
    maxLength <- max(length(allChildren), length(allParents))
    allChildren <- c(allChildren, rep("", maxLength - length(allChildren)))
    allParents <- c(allParents, rep("", maxLength - length(allParents)))
    df <- data.frame(Parents = allParents,
                     Children = allChildren,  
                     stringsAsFactors = FALSE)
    df
  })

  output$SummaryGraph <- renderForceNetwork({
    if(is.null(P_selection())) return(NULL)
    
    if(is.null(SummGraph())){
      Links_sum <- data.frame(source = 0, target = 0, value = 1)
      Nodes_sum <- data.frame(name = prot_names_short[P_selection()], group = "Selected", size = 1)
    }else{
      Links_sum <- SummGraph()$Links_sum
      Nodes_sum <- SummGraph()$Nodes_sum
    }
    fN <- forceNetwork(Links = Links_sum, Nodes = Nodes_sum,
                 Source = "source", Target = "target",
                 Value = "value", NodeID = "name",
                 Group = "group", opacity = 0.99, 
                 arrows = T, zoom = T, charge = -20,
                 opacityNoHover = TRUE, legend = T,
                 colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
    fN
  })

  output$TemporalGraph <- renderForceNetwork({
    if(is.null(P_selection())) return(NULL)
    
    fN <- forceNetwork(Links = TempGraph()$Links_temp, Nodes = TempGraph()$Nodes_temp,
                       Source = "source", Target = "target",
                       Value = "value", NodeID = "name",
                       Group = "group", opacity = 0.99,# Nodesize = 3,
                       arrows = T, zoom = T, legend=T, charge = -15,
                       opacityNoHover = TRUE,
                       linkColour = TempGraph()$direction,
                       colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
    # Inject custom D3 JavaScript with Capture-Phase Safety Windows
    fN <- htmlwidgets::onRender(fN, '
      function(el, x) {
        var nodes = d3.select(el).selectAll(".node").data();
        
        // 1. Create the UI panel and buttons
        var panel = document.createElement("div");
        panel.style.position = "absolute";
        panel.style.bottom = "20px";
        panel.style.left = "50%";
        panel.style.transform = "translateX(-50%)";
        panel.style.display = "flex";
        panel.style.gap = "15px";
        panel.style.zIndex = "1000";
        el.appendChild(panel);
        
        function styleButton(btn, text, bgColor) {
          btn.innerHTML = text;
          btn.style.padding = "10px 18px";
          btn.style.fontSize = "14px";
          btn.style.fontWeight = "bold";
          btn.style.cursor = "pointer";
          btn.style.borderRadius = "6px";
          btn.style.border = "1px solid #ccc";
          btn.style.backgroundColor = bgColor;
          btn.style.boxShadow = "0px 2px 4px rgba(0,0,0,0.1)";
          
          // Stop accidental background graph dragging
          btn.addEventListener("mousedown", function(e) { e.stopPropagation(); });
          btn.addEventListener("touchstart", function(e) { e.stopPropagation(); });
          
          panel.appendChild(btn);
        }
        
        var btnSep = document.createElement("button");
        styleButton(btnSep, "Separate by time-point", "#fff");
        
        var btnRelax = document.createElement("button");
        styleButton(btnRelax, "Relax", "#fff");
        
        var btnFree = document.createElement("button");
        styleButton(btnFree, "Free Up", "#e0f7fa"); 
        
        // 2. State Tracking and Memory
        var isSeparating = false;
        var isRelaxing = false;
        var initialPositions = {};
        var memorySaved = false;
        
        var dragTarget = null;
        var dragX = 0;
        var dragY = 0;
        
        function saveMemory() {
          if (!memorySaved) {
            nodes.forEach(function(d) {
              initialPositions[d.index] = { x: d.x, y: d.y };
            });
            memorySaved = true;
          }
        }
        
        // 3. Persistent Wake Engine using D3 Drag State
        function sendDragStart() {
          try {
            dragTarget = el.querySelector(".node circle") || el.querySelector(".node");
            if (dragTarget) {
              var rect = dragTarget.getBoundingClientRect();
              dragX = rect.left + rect.width / 2;
              dragY = rect.top + rect.height / 2;
              
              dragTarget.dispatchEvent(new MouseEvent("mousedown", {
                bubbles: true, cancelable: true, view: window, button: 0, clientX: dragX, clientY: dragY
              }));
            }
          } catch(e) {}
        }
        
        function sendDragEnd() {
          try {
            if (dragTarget) {
              window.dispatchEvent(new MouseEvent("mouseup", {
                bubbles: true, cancelable: true, view: window, button: 0, clientX: dragX, clientY: dragY
              }));
              dragTarget = null;
            }
          } catch(e) {}
        }
        
        // 4. Interaction Handlers
        function startInteraction(type) {
          saveMemory();
          if (type === "sep") isSeparating = true;
          if (type === "relax") isRelaxing = true;
          
          // Pre-lock nodes to current position so they do not jump when beam engages
          nodes.forEach(function(d) {
            if (d.group === "6h" || d.group === "48h") {
              d.fx = d.x;
              d.fy = d.y;
            }
          });
          
          sendDragStart();
        }
        
        function stopInteraction() {
          if (!isSeparating && !isRelaxing) return;
          isSeparating = false;
          isRelaxing = false;
          
          sendDragEnd();
          
          // Nodes remain cleanly pinned exactly where you let go
          nodes.forEach(function(d) {
            if (d.group === "6h" || d.group === "48h") {
              if (d.fx !== null) d.x = d.fx;
              if (d.fy !== null) d.y = d.fy;
            }
          });
        }
        
        // 5. Global Capture Phase Listeners (Master Kill Switch)
        btnSep.addEventListener("mousedown", function() { startInteraction("sep"); });
        btnSep.addEventListener("touchstart", function(e) { e.preventDefault(); startInteraction("sep"); });
        
        btnRelax.addEventListener("mousedown", function() { startInteraction("relax"); });
        btnRelax.addEventListener("touchstart", function(e) { e.preventDefault(); startInteraction("relax"); });
        
        // The true argument forces capture phase execution to intercept swallowed events
        window.addEventListener("mouseup", stopInteraction, true);
        window.addEventListener("touchend", stopInteraction, true);
        window.addEventListener("blur", stopInteraction, true);
        
        // 6. Free Up Button
        btnFree.addEventListener("click", function() {
          isSeparating = false;
          isRelaxing = false;
          
          nodes.forEach(function(d) { d.fx = null; d.fy = null; });
          
          sendDragStart();
          setTimeout(function() {
            sendDragEnd();
            nodes.forEach(function(d) { d.fx = null; d.fy = null; });
          }, 300);
        });
        
        // 7. Custom Physics Loop (Tractor Beam Animation)
        d3.timer(function() {
          if (!isSeparating && !isRelaxing) return;
          
          var width = el.clientWidth;
          var speed = 0.05; 
          
          // Micro-pulse the drag to forcefully keep the main D3 simulation awake
          if (dragTarget) {
            dragX += 0.1;
            window.dispatchEvent(new MouseEvent("mousemove", {
              bubbles: true, cancelable: true, view: window, button: 0, clientX: dragX, clientY: dragY
            }));
          }
          
          nodes.forEach(function(d) {
            if (d.group === "6h" || d.group === "48h") {
              
              // Failsafe initialization
              if (d.fx == null) d.fx = d.x;
              if (d.fy == null) d.fy = d.y;
              
              if (isSeparating) {
                var targetX = (d.group === "6h") ? width * 0.1 : width * 0.9;
                // Animate absolute locks instead of weak velocity
                d.fx += (targetX - d.fx) * speed;
                
              } else if (isRelaxing) {
                var origX = initialPositions[d.index].x;
                var origY = initialPositions[d.index].y;
                
                d.fx += (origX - d.fx) * speed;
                d.fy += (origY - d.fy) * speed;
              }
            }
          });
        });
      }
    ')
    
    fN
  })

  output$downloadPvalsDrug <- downloadHandler(
    filename = function() {
      paste('DrugEffects-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      if(is.null(P_selection())) {
        df_res <- data.frame(Drug_or_Interaction = "No protein selected!", PValue = NA, 
                             stringsAsFactors = FALSE)
      }else{
        df_res <- data.frame(Drug_or_Interaction = names(pvec()), PValue = pvec(), 
                            stringsAsFactors = FALSE)
      }
      write.csv(df_res, con, row.names = FALSE)
  })

  output$downloadPvalsProtein <- downloadHandler(
    filename = function() {
      paste('ProteinEffects-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      Links_all_res <- Links_all()
      
      # Safety Check 1: Is there any data?
      if(is.null(P_selection()) || is.null(Links_all_res) || length(Links_all_res) == 0) {
        df_res <- data.frame(source = "No data or no protein selected!", 
                             target = NA, PValue = NA, stringsAsFactors = FALSE)
      } else {
        # Safety Check 2: Use seq_along to avoid the 1:0 trap
        df_res <- do.call(rbind, lapply(seq_along(Links_all_res), function(i){
          links <- Links_all_res[[i]]
          
          # Skip if this specific transition is empty
          if(is.null(links) || nrow(links) == 0) return(NULL)
          
          # Ensure we don't go out of bounds of expTimes
          if (i + 1 > length(expTimes)) {
            transition_label <- paste0(expTimes[i], "h to unknown")
          } else {
            transition_label <- paste0(expTimes[i], "h to ", expTimes[i+1], "h")
          }
          
          data.frame(InteractionType = transition_label,
                     Source = prot_names_short[links[, "source"]], 
                     Target = prot_names_short[links[, "target"]], 
                     PValue = links[, "pvalue"], 
                     stringsAsFactors = FALSE)
        }))
      }
      
      # Handle cases where all transitions were NULL/empty
      if(is.null(df_res)) {
        df_res <- data.frame(Status = "No significant interactions found")
      }
      
      write.csv(df_res, con, row.names = FALSE)
    }
  )

  output$downloadSummary <- downloadHandler(
    filename = function() {
      paste('ProteinNetworkSummary-', Sys.Date(), '.html', sep='')
    },
    content = function(con) {
      if(is.null(P_selection())) return(NULL)
      
      if(is.null(SummGraph())){
        Links_sum <- data.frame(source = 0, target = 0, value = 1)
        Nodes_sum <- data.frame(name = prot_names_short[P_selection()], group = "Selected", size = 1)
      }else{
        Links_sum <- SummGraph()$Links_sum
        Nodes_sum <- SummGraph()$Nodes_sum
      }
      fN <- forceNetwork(Links = Links_sum, Nodes = Nodes_sum,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   Group = "group", opacity = 0.99, 
                   arrows = T, zoom = T, charge = -20,
                   opacityNoHover = TRUE, legend = T,
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
      saveNetwork(fN, file = con)
  })

  output$downloadTemporal <- downloadHandler(
    filename = function() {
      paste('ProteinNetworkTemporal-', Sys.Date(), '.html', sep='')
    },
    content = function(con) {
      if(is.null(P_selection())) return(NULL)
      fN <- forceNetwork(Links = TempGraph()$Links_temp, Nodes = TempGraph()$Nodes_temp,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   Group = "group", opacity = 0.99,# Nodesize = 3,
                   arrows = T, zoom = T, legend=T, charge = -15,
                   opacityNoHover = TRUE,
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
      saveNetwork(fN, file = con)
  })
  
  # Automatically select the treatment with the lowest p-value whenever pvec updates
  observeEvent(pvec(), {
    req(pvec()) # Ensure pvec is calculated and not NULL
    
    # 1. Find the actual text name of the drug combination with the minimum p-value
    # Using names() is much safer than positional indexing (like treatNames[which.min(...)])
    best_treatment <- treatNames[which.min(pvec())]
    
    # 2. Update the UI drop-down menu on the fly
    updateSelectInput(
      session = getDefaultReactiveDomain(), 
      inputId = "selTreat", 
      selected = best_treatment
    )
  })
}


shinyApp(ui, server)
