
## Need to add doc here.

custom_moves <- function(...) {
    
    move_functions <- list(...)
    
    if (length(move_functions) == 1L && is.list(move_functions[[1]])) {
        move_functions <- move_functions[[1]]
    }
    
   
    defaults <- list(mu = cpp_move_mu,
                     pi = cpp_move_pi,
                     alpha = cpp_move_alpha,
                     swap_cases = cpp_move_swap_cases,
                     t_inf = cpp_move_t_inf,
                     kappa = cpp_move_kappa
                     )

     
    moves <-  modify_defaults(defaults, move_functions, FALSE)
    moves_names <- names(moves)

    
    
    ## check all moves are functions

    function_or_null <- function(x) {
        is.null(x) || is.function(x)
    }
    
    is_ok <- vapply(moves, function_or_null, logical(1))
    
    if (!all(is_ok)) {
        culprits <- moves_names[!is_ok]
        msg <- paste0("The following moves are not functions: ",
                      paste(culprits, collapse = ", "))
        stop(msg)
    }

    
    ## check they all have a 'param' argument

    param_is_arg <- function(x) {
        if(is.function(x)) {
            return ("param" %in% methods::formalArgs(x))
        }
        
        return(TRUE)
    }
    
    param_ok <- vapply(moves, param_is_arg, logical(1))

    if (!all(param_ok)) {
        culprits <- moves_names[!param_ok]
        msg <- paste0("The following moves dont' have a 'param' argument: ",
                      paste(culprits, collapse = ", "))
        stop(msg)
    }
    

    return(moves)

}
