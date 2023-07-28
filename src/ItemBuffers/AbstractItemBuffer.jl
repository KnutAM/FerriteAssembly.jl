abstract type AbstractItemBuffer end

# Acess functions
function get_ae end
function get_aeold end
function get_re end
function get_Ke end
function get_material end
function get_values end
function get_old_state end
function get_new_state end 
function get_time_increment end 
function get_user_data end 
function get_user_cache end 

#Ferrite.celldofs(::AbstractItemBuffer)
#Ferrite.getcoordinates(::AbstractItemBuffer)
#Ferrite.cellid(::AbstractItemBuffer) 
#Ferrite.dof_range(::AbstractItemBuffer, ::Symbol)

# Set functions  
function set_time_increment! end

# Common parts for taskLocals interface
function scatter!(task::AbstractItemBuffer, base::AbstractItemBuffer)
    set_time_increment!(task, get_time_increment(base))
end
gather!(::AbstractItemBuffer, ::AbstractItemBuffer) = nothing