#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyalert)
library(shinyjs)
library(combinat)
library(pracma)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#set.seed(1) 

# Génére une répartition

divisors <- function(x){
  #  Vector of numbers to test against
  y <- seq_len( x)
  #  Modulo division. If remainder is 0 that number is a divisor of x so return it
  res = y[ x%%y == 0 ]
  cbind(res,rev(res))
}

write_divisors <- function(x){
  l<-list()
  for (i in 1:nrow(x)) {
    l = append(l,  paste(x[i,1],"x", x[i,2]))
  }
  return(l)
}

find_selected_div <- function(s,l){
  if (s %in% l){
    for (i in 1:length(l)) {
      if (s %in% l[i]){
        return(i)
      }
    }
  }
}

rc_ind <- function(M, ind) c(row(M)[ind], col(M)[ind] )

ind_from_block_id <- function(id,block_mat=block_mat,block_dim=block_dim){
  ind_tmp = rc_ind(block_mat,id)
  print(ind_tmp)
  row = (block_dim[1]*(ind_tmp[1]-1)+1):((block_dim[1]*(ind_tmp[1]-1)+1)+block_dim[1]-1)
  col =  (block_dim[2]*(ind_tmp[2]-1)+1):((block_dim[2]*(ind_tmp[2]-1)+1)+block_dim[2]-1)
  list(row=row, col=col)
}

create_constraint_list <- function(M,ind_sub,constr_window){
  fail=0
  max_row = nrow(M)
  max_col = ncol(M)
  min_row = 1
  min_col = 1
  constraint_list = list()
  #print(M)
  #print(ind_sub)
  if(isempty(constr_window)){
    constraint_list <- vector("list", length(ind_sub$row)*length(ind_sub$col))
  }
  else{
    
    for (i in ind_sub$row) {
      
      for (j in ind_sub$col) {
        
        #   print("i,j")
        #  print(i)
        #  print(j)
        constr_window_list = c()
        for (k in 1:nrow(constr_window)){
          #    print("k")
          #   print(k)
          current_window_row = i+constr_window[k,1];
          current_window_col = j+constr_window[k,2]; 
          # print(current_window_row)
          # print(current_window_col)
          
          if ((current_window_row>=min_row) & (current_window_row<=max_row) & (current_window_col>=min_col) & (current_window_col<=max_col) ){
            # si un point est déjà présent avec une valeur positive, on l'ajoute à la liste des contraintes temporaire
            #  print("tutu")
            #  print(c(current_window_row,current_window_col))
            #  print( M[current_window_row,current_window_col])
            if (M[current_window_row,current_window_col]>0){
              constr_window_list=append(constr_window_list,M[current_window_row,current_window_col])
            }
          }
          
        }
        if (isempty(constr_window_list)){
          constr_window_list=append(constr_window_list,-1)
        }
        constraint_list=append(constraint_list,list(constr_window_list))
      }
    }
  }
  return(constraint_list)
  
}

block_generation <- function(constraint_list){
  # on trie la liste pour avoir le point avec le plus de contraintes en premier et faciliter la génération
  if ( !isempty(unlist(constraint_list))){
    sorted_ind =  order(sapply(constraint_list,'[[',1),decreasing = TRUE)
    modalites = 1:length(sorted_ind)
    available_mod = rep(1,length(sorted_ind))
    lin_block = rep(0,length(sorted_ind))
    fail_flag=0
    
    for (i in 1:length(sorted_ind)) {
      
      current_ind = sorted_ind[i]
      
      tmp_av_mod = available_mod
      if (length(constraint_list[[current_ind]]==1) & (constraint_list[[current_ind]][1]==-1)){}
      else{
        tmp_av_mod[constraint_list[[current_ind]]]=0
      }
      
      if(sum(modalites[tmp_av_mod==1])==0){
        fail_flag=1 
      }
      else{
        if (length(modalites[tmp_av_mod==1])==1){
          lin_block[current_ind] = modalites[tmp_av_mod==1]
          available_mod[ modalites[tmp_av_mod==1]]=0
        }
        else{
          lin_block[current_ind] = sample(modalites[tmp_av_mod==1],1)
          available_mod[ lin_block[current_ind]]=0
        }
      }
      
    }
  }
  else
  {
    # cas sans contrainte, on génére juste aléatoirement un bloc
    lin_block = sample(1:length(constraint_list),length(constraint_list))
    fail_flag=0
  }
  
  list(lin_block=lin_block, fail_flag=fail_flag)
}


# numericInput("n_modalites", "Number of crop modalities", value = 6, min = 1, max = 8, step = 1),
# selectInput("block_shape", "Block shape (rows x columns)", c("")),
# numericInput("n_rep", "Number of replicates", value = 2, min = 1, max = 20, step = 1),
# selectInput("plot_shape", "Plot shape (rows x columns)", c("")),
# selectInput("constraint", "Constraint", c("None", "Orthogonal","Ortho+Diag","Orthogonal 2","Ortho2 + Diag2")),

calculate_repartition <- function(n_modalites,block_shape,n_rep,plot_shape,constraint) {
  
  # récupération des dimensions de bloc et de parcelle
  possible_block_dim  = divisors(n_modalites);
  list_block_dim = write_divisors(possible_block_dim)
  possible_plot_dim  = divisors(n_rep); 
  list_plot_dim = write_divisors(possible_plot_dim)
  
  ind_block_div = find_selected_div(block_shape,list_block_dim)
  ind_plot_div = find_selected_div(plot_shape,list_plot_dim)
  
  
  block_dim = possible_block_dim[ind_block_div,]
  parc_block_dim= possible_plot_dim[ind_plot_div,]
  
  # Contraintes...
  empty_window = {}
  ortho_window = cbind(c(-1,1,0,0),c(0,0,-1,1))
  ortho_diag_window = cbind(c(-1:1,-1,1,-1:1),c(-1,-1,-1,0,0,1,1,1))
  ortho_window_2 = rbind(ortho_diag_window, cbind(c(-2,2,0,0),c(0,0,-2,2)))
  ortho_diag_window_2 = cbind(c(-2:2,-2:2,-2:2,-2:2,-2:2),c(-2,-2,-2,-2,-2,-1,-1,-1,-1,-1,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2))
  neighb_windows = list(empty_window,ortho_window, ortho_diag_window,ortho_window_2,ortho_diag_window_2)
  constr_list =  c("None", "Orthogonal","Ortho+Diag","Orthogonal 2","Ortho2 + Diag2")
  
  ind_window =find_selected_div( constraint,constr_list)
  
  
  constr_window = neighb_windows[ind_window][[1]]
  
  
  # Matrice finale
  M <- matrix(0, nrow=parc_block_dim[1]*block_dim[1], ncol=parc_block_dim[2]*block_dim[2])
  
  # Matrice des blocs
  block_id = 1:(parc_block_dim[1]*parc_block_dim[2])
  block_mat  = matrix(0, nrow=parc_block_dim[1], ncol=parc_block_dim[2])
  block_mat[]=block_id
  
  # création d'une matrice pour l'affichage 
  ind = 1;
  id_mat<- matrix(0, nrow=parc_block_dim[1]*block_dim[1], ncol=parc_block_dim[2]*block_dim[2])
  for (i in 1:nrow(block_mat)) {
    for (j in 1:ncol(block_mat)) {
      ind_current_block = ind_from_block_id(ind,block_mat,block_dim)
      id_mat[ind_current_block$row,ind_current_block$col] = ind
      ind = ind +1
    }  
  }
  
  
  # Calculer toutes les permutations possibles des nombres 1 à n_modalites
  #permutations <- permn(1:n_modalites)
  
  # le premier bloc n'a pas de contraintes
  ind_sub=ind_from_block_id(1,block_mat,block_dim )
  M[ind_sub$row,ind_sub$col]<-1:n_modalites
  
  
  
  # nombre max d'iterations
  iter_max = 10000
  internal_iter_max = 100
  
  # compteur d'iterations et de blocs
  n_iter <- 1
  n_block = 2
  
  # critere d'arret
  end_alg =0
  fail_flag = 0
  while (n_iter < iter_max & end_alg==0){
    # on créé les indices du nouveau bloc
    ind_sub=ind_from_block_id(n_block,block_mat,block_dim )
    
    # on créé la liste des contraintes associées
    constraint_list = create_constraint_list(M,ind_sub,constr_window)
    
    #on génère le nouveau bloc
    block_found = 0
    n_internal_iter = 0
    while (n_internal_iter < internal_iter_max & block_found==0){
      block_gen= block_generation(constraint_list)
      
      # si ça a trouvé un bloc respectant les contraintes, on arrete, sinon on recommence
      if (block_gen$fail_flag ==0){
        block_found = 1
        lin_block=block_gen$lin_block 
        
        ind = 1
        for (i in ind_sub$row) {
          for (j in ind_sub$col) {
            M[i,j]<-lin_block[ind]
            ind = ind +1;
          }
        }
        # peut etre plus propre mais ne fonctionne pas
        #M[ind_sub$row,ind_sub$col]<-lin_block
      }
      n_internal_iter = n_internal_iter +1
      n_iter = n_iter + 1
    }
    
    if (block_found==0 | (block_found==1 & n_block== n_rep)){
      #si la contrainte ne peut pas être respectée, on essaye d'enlever le bloc d'avant
      if (block_found==0  & n_internal_iter == internal_iter_max & n_block>2){
        if (n_block> parc_block_dim[2]){
          #pour éviter d'être complètement bloqué à cause d'un bloc antérieur, on enlève jusqu'au n_block- parc_block_dim(2)
          
          for (k in 1:parc_block_dim[2])
          {
            ind_sub=ind_from_block_id(n_block-1,block_mat,block_dim )
            
            for (i in ind_sub$row) {
              for (j in ind_sub$col) {
                M[i,j]<-0
              }
            }
            n_block = n_block-1
          }
          
        }
        #on essaye d'en enlever un seul dans la plupart des cas
        else{
          ind_sub=ind_from_block_id(n_block-1,block_mat,block_dim )
          
          for (i in ind_sub$row) {
            for (j in ind_sub$col) {
              M[i,j]<-0
            }
          }
          n_block = n_block-1
        }
        
        
        
      }
      else{
        end_alg=1
      }
      
    }
    else{
      n_block = n_block+1
    }
    
  }
  
  if (n_iter==iter_max |block_found==0)
  {fail_flag = 1}
  
  
  list(M=M, id_mat=id_mat,n_iter = n_iter,fail_flag = fail_flag, block_mat = block_mat, block_dim = block_dim)
}













# Define server logic
function(input, output, session) {

  observe({
    n_modalites <- input$n_modalites
    
    if (n_modalites >20){ updateNumericInput(session,"n_modalites",value=20)}
    if (n_modalites <1){ updateNumericInput(session,"n_modalites",value=1)}
    
    # Liste des dimensions possibles de blocs
    possible_block_dim  = divisors(n_modalites);
    
    # Can also set the label and select items
    updateSelectInput(session, "str_block_shape",
                      choices = write_divisors(possible_block_dim),
                      selected = 1)
  })
  
  observe({
    n_rep <- input$n_rep
    
    if (n_rep >20){ updateNumericInput(session,"n_rep",value=20)}
    if (n_rep <1){ updateNumericInput(session,"n_rep",value=1)}
    
    # Liste des dimensions possibles de blocs
    possible_parc_dim  = divisors(n_rep);
    
    # Can also set the label and select items
    updateSelectInput(session, "str_plot_shape",
                      choices = write_divisors(possible_parc_dim),
                      selected = 1)
  }) 
  
  
  output_change <- eventReactive(input$calculate, {
    print(input$block_shape)
    shinyjs::disable("calculate")
    
    # repartition.output = calculate_repartition(input$n_c, input$n_blocs_col, input$n_modalites, input$diag_cont,as.numeric(input$iter_max) )})
    withProgress(message = 'Doing Simulations...', value = 1, {
      repartition.output = calculate_repartition(input$n_modalites,input$str_block_shape,input$n_rep,input$str_plot_shape,input$constraint )
      shinyjs::enable("calculate")
   
      list(M=repartition.output$M, id_mat=repartition.output$id_mat, n_iter=repartition.output$n_iter,fail_flag =repartition.output$fail_flag , block_mat=repartition.output$block_mat,block_dim =repartition.output$block_dim )
 
       })
  })
  
  
  
  # Affichage 
  output$repartition_plot <- renderPlot({
    M<-output_change()$M
    id_mat<-output_change()$id_mat
    n_iter<-output_change()$n_iter
    fail_flag<-output_change()$fail_flag  
    block_mat = output_change()$block_mat
    block_dim = output_change()$block_dim
    if(fail_flag ==1){
      shinyalert("Oops!", "Not found, try again, increase the number of iterations or release the constraints.", type = "error")
    }else{
 
      M=t(M)
      id_mat = t(id_mat)
      image(1:nrow(M),1:ncol(M),(M),xaxt = "n",xlab="Columns",ylab="Rows",col=hcl.colors(input$n_modalites, "terrain",alpha=0.8))
      text(x = row((M)),
           y = col((M)),
           # label = sprintf("%d: %d", (id_mat),(M)))
           label = sprintf("%d", (M)))
      axis(1, at = 1:nrow(M))

      # set the colour palette
      cols <- hcl.colors(input$n_rep*1.7,"Plasma")
    
      for (i in 1:input$n_rep) {

        res = ind_from_block_id(i,block_mat,block_dim)
        par(new=TRUE)
        rect(xleft = min(res$col)-0.48, xright =max(res$col)+0.48, ybottom = min(res$row)-0.48, ytop = max(res$row)+0.48, lwd =7,col=cols[i],density=0.)
        text(mean(res$col), max(res$row)+0.2, labels = paste("Bloc", i), cex = 1.5,col=cols[i])
      }
    }
  }
  )
  
  # Exportation au format CSV
  output$downloadCSV <- downloadHandler(
    filename = function() {
      paste("donnees_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(output_change()$id_mat*100+output_change()$M, file)
      
    }
  )
}
