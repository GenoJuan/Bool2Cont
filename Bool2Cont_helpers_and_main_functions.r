library(ggplot2)
library(BoolNet)
library(reshape)

#####################################################################################################################
#####################################################################################################################

node.rm = function(net, node2rm = "random", exclude = c(), verbose = TRUE){
    
    ############################ Select a node to remove ##########################################
    
    if(node2rm == "random"){
        numNodes = length(net$genes)
        nodes = seq(1, numNodes)
        node.i = sample(nodes, 1)
    }
    else if(class(node2rm) == "numeric"){
        if(node2rm < 1 || node2rm > length(net$genes)){
            stop('Supplied numeric \"node2rm\". node2rm must be in the range from 1 to net$genes length')
        } else{
            numNodes = length(net$genes)
            node.i = node2rm
        }
    } 
    else if(class(node2rm) == "character"){
        numNodes = length(net$genes)
        node.i = which(net$genes == node2rm)
        if(length(node.i) == 0){
            stop(paste("Error: No gene named", node2rm, "in net object."))
        }
    }
        else {
        stop('Error type: \"node2rm\". node2rm must be either "random" or numeric.')
    }
        
    ##############################################################################################
    #################################### Generate a network copy #################################
    
    net_tmp = net
    
    ##############################################################################################
    #################### Retrieve the Boolean function for the selected node #####################
        
    BF.i = net$interactions[[node.i]]$expression
        
    ##############################################################################################
    ######### In every BF of nodes regulaed by i, sustitute the "node i" term by the BF i ########
    
    for(node in 1:numNodes){
        if(length(which(net$interactions[[node]]$input == node.i)) == 1){
          net_tmp$interactions[[node]]$expression = gsub(net_tmp$genes[node.i], paste("(", BF.i, ")"), net$interactions[[node]]$expression, perl = TRUE)
        }
    }
        
    ##############################################################################################
    #################################### Write network in a tmp file #############################
    
    sink("net_tmp.net")
    cat("targets, factors\n")
    
    for(node in 1:numNodes){
        if(node != node.i){
            exp = paste(net_tmp$genes[node], ",", net_tmp$interactions[[node]]$expression, "\n")
            cat(exp)
          }
    }
    
    sink()
        
    ##############################################################################################
    ######################################### Process network file ###############################
    
    #return(BF.i)
    net_tmp = tryCatch({
        net_tmp = loadNetwork("net_tmp.net")
        saveNetwork(net_tmp, "net_temp.net", generateDNFs = "short")
        net_tmp = loadNetwork("net_temp.net")
        return(net_tmp)
        }, 
        warning = function(x){net})
    system("rm net_tmp.net")
    return(net_tmp)
}
  
 
 
  
  



###################################################################################################
###################################################################################################