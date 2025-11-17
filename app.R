## app.R ##
library(shinydashboard)
library(plotly)
library(networkD3)
library(grid)
library(readxl)
library(shinyWidgets)

#choises of timepoints for drug effects
t_choice <- c("6h", "24h", "48h")

# replace drug ids by names
load('data/drugLookup.RData')
replace_drug_ids <- function(x) {
  ids <- strsplit(x, "_|`")[[1]]
  ids <- ids[grepl("#", ids)]
  names <- sapply(ids, function(id) drug_lookup[[id]])
  paste(names, collapse = ":")
}

# load ordering of drugs (sort drugs with experiments together)
load('data/order.RData')
nDrugs <- length(drugOrder)
drugOrder <- sapply(drugOrder, replace_drug_ids)

ui <- dashboardPage(
  dashboardHeader(title = "Drug Perturbations"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Settings", tabName = "Settings", icon = icon("cog")),
      menuItem("Drug Effects", tabName = "DrugEffects", icon = icon("dashboard")),
      menuItem("Protein Network", tabName = "ProteinNetwork", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "Settings", 
              h2("Settings"),
              fluidRow(
                box(title = "Select Proteins of Interest", status = "primary", solidHeader = TRUE,
                    pickerInput("protSet", "Select Protein Set", 
                                choices = unname(prot_names_short), 
                                multiple = TRUE, 
                                options = list("live-search"=TRUE)),
                    actionButton("preSelected", "pre Selected Set"),
                    actionButton("clear", "Clear Selection"), 
                    fileInput("file", "Upload .txt File with Protein Names (one per line)", accept = c(".txt")),
                    h3("p-value Adjustments"),
                    numericInput("alpha", "Significance Level", value = 0.05, min = 0, max = 1, step = 0.0001),
                    selectInput("corectionDrug", "Correction Method for Drug Effects", choices = p.adjust.methods, selected = 'holm'),
                    selectInput("corectionProtein", "Correction Method for Protein Effects", choices = p.adjust.methods, selected = 'BH'),
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
              fluidRow(
                box(title = "Drug Effect Heatmap", status = "primary", solidHeader = TRUE,
                    plotlyOutput("plotDrugEffects", height = 800), width = 6),
                box(title = "Significant Single Drugs", status = "primary", solidHeader = TRUE,
                    tableOutput("singleDrugs"), width = 3),
                box(title = "Significant Drug Interactions", status = "primary", solidHeader = TRUE,
                    tableOutput("interactions"), width = 3)
              )
      ),
      tabItem(tabName = "ProteinNetwork", 
              h2("Protein Network"),
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
  
  # load data for Drug Effects
  load("data/DrugEffects.RData")
  load("data/proteinSelection.RData")
  load("data/proteinNetworkPval.RData")


  observeEvent(input$preSelected, {
    updatePickerInput(session = getDefaultReactiveDomain(), inputId = "protSet", 
                      selected = prot_names_short[prot_names_short %in% most_imp])
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
  P_selection <- reactive({
    sel <- which(prot_names_short %in% input$protSet)
    if(length(sel) == 0) return(NULL)
    sel
  })

  pvec <- reactive({
    if(is.null(input$t)) return(NULL)
    t_selection <- t_choice %in% input$t
    
    #adjust p values of drug effects on selected proteins
    selPvecs <- array(allPvecs[, , P_selection()], dim = c(dim(allPvecs)[1:2], length(P_selection())))
    selPvecs <- array(p.adjust(selPvecs, method = input$corectionDrug), dim = dim(selPvecs))
    
    # collect min p value of drug effect over proteins and time points
    pvec <- apply(selPvecs, 1, function(p) min(p[t_selection, ]))
    names(pvec) <- sapply(treatment, replace_drug_ids)
    pvec
  })
  
  Links_all <- reactive({
    if(is.null(P_selection())) return(NULL)
    
    # select relevant p-values
    Pval_sel <- lapply(Pval_all, function(links){
      rel.Links <- links$source %in% P_selection() | links$target %in% P_selection()
      links[rel.Links, ]
    })
    
    # p-value correction
    Tpval <- sapply(Pval_sel, function(links) links$pvalue)
    Tpval.corr <- matrix(p.adjust(Tpval, method = input$corectionProtein), ncol = 2)
    Pval_sel[[1]][, "pvalue"] <- Tpval.corr[, 1]
    Pval_sel[[2]][, "pvalue"] <- Tpval.corr[, 2]
    
    # convert p-values to links
    Links_all <- lapply(Pval_sel, function(links) links[links$pvalue < input$alpha, ])
    
    if(sum(sapply(Links_all, nrow)) == 0) return(NULL)
    Links_all
  })

  SummGraph <- reactive({
    if(is.null(Links_all())) return(NULL)
    Links_sum <- do.call(rbind, Links_all())
    Links_sum$source <- Links_sum$source - 1
    Links_sum$target <- Links_sum$target - 1
    Links_sum$value <- 1
  
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
    #if(is.null(Links_all())) return(NULL)
    expTimes <- c(6, 24, 48)
    rel6 <- sort(unique(c(P_selection(), Links_all()[[1]][, "source"])))
    rel24 <- sort(unique(c(P_selection(), Links_all()[[1]][, "target"], Links_all()[[2]][, "source"])))
    rel48 <- sort(unique(c(P_selection(), Links_all()[[2]][, "target"])))
    rel <- list(rel6, rel24, rel48)
    lenRel <- c(0, length(rel6), length(rel24), length(rel48))
    
    nodenames <- c(paste(prot_names_short[rel6], expTimes[1], sep = '_'), 
                   paste(prot_names_short[rel24], expTimes[2], sep = '_'), 
                   paste(prot_names_short[rel48], expTimes[3], sep = '_'))
    nodegroups <- rep(paste0(expTimes, "h"), times = c(length(rel6), length(rel24), length(rel48)))
    
    if(is.null(Links_all())){
      Links_temp <- data.frame(source = 0, target = 0, value = 1)
    }else{
      Links_temp <- lapply(1:2, function(t){
        links <- Links_all()[[t]]
        for(i in 1:length(rel[[t]])){
          links[links[, 1] == rel[[t]][i], 1] <- i - 1 + sum(lenRel[1:t])
        }
        for(i in 1:length(rel[[t+1]])){
          links[links[, 2] == rel[[t+1]][i], 2] <- i - 1 + sum(lenRel[1:(t+1)])
        }
        links
      })
      Links_temp <- do.call(rbind, Links_temp)
      Links_temp$value <- 1
    }
    
    Nodes_temp <- data.frame(name = nodenames, group = nodegroups, size = 0.3)
    Nodes_temp$radius <- as.numeric(c(rel6, rel24, rel48))
    list(Links_temp = Links_temp, Nodes_temp = Nodes_temp)
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
      Children <- sum(SummGraph()$Links_sum$source == i)#-1
      Parents <- sum(SummGraph()$Links_sum$target == i)#-1
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
      Children <- SummGraph()$Nodes_sum$name[SummGraph()$Links_sum$target[SummGraph()$Links_sum$source == i & SummGraph()$Links_sum$target != i] + 1]
      Parents <- SummGraph()$Nodes_sum$name[SummGraph()$Links_sum$source[SummGraph()$Links_sum$target == i & SummGraph()$Links_sum$source != i] + 1]
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
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
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
      if(is.null(P_selection())) {
        df_res <- data.frame(source = "No protein selected!", target = NA, PValue = NA, 
                             stringsAsFactors = FALSE)
      }else{
        Links_all_res <- Links_all()
        df_res <- do.call(rbind, lapply(1:length(Links_all_res), function(i){
          links <- Links_all_res[[i]]
          if(nrow(links) == 0) return(NULL)
        df <- data.frame(InteractionType = ifelse(i == 1, "6h to 24h", "24h to 48h"),
                         Source = prot_names_short[links$source], 
                         Target = prot_names_short[links$target], 
                         PValue = links$pvalue, 
                         stringsAsFactors = FALSE)
          df
        }))
      }
      write.csv(df_res, con, row.names = FALSE)
  })

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
}


shinyApp(ui, server)