## Memory handling for memory

#
struct CodeMem
    slots::Vector{Symbol}
    slot_names::Vector{String}
    special_names::Dict{Symbol,String}
end

# Creates a memory object
# k is the number of slots
# slot_name_fun(i) should return a string
#    which is the name of the node
function CodeMem(k, slot_name_fun)
    slots = map(i -> :Free, 1:k) # Special keyword
    #slot_names=Vector{String}(undef,k);
    slot_names = slot_name_fun.(1:k)
    return CodeMem(slots, slot_names, Dict{Symbol,String}())
end

function get_slot_number(mem::CodeMem, node::Symbol)
    return findfirst(mem.slots .== node)
end

function set_slot_number!(mem, slot_number, node)
    return mem.slots[slot_number] = node
end

function get_slot_name(mem::CodeMem, node::Symbol)
    if (haskey(mem.special_names, node))
        return mem.special_names[node]
    else
        if (node in mem.slots)
            return mem.slot_names[findfirst(mem.slots .== node)]
        else
            @show mem.slots
            error("Unable to find memory slot for $node")
        end
    end
end

function is_allocated(mem::CodeMem, node::Symbol)
    return node in mem.slots
end

# Returns a reference to a memory slot
# which is either a new slot or
# one of the elements in the Vector allow_inplace
#
# Returns a Tuple with a reference to the slot
# and the name
#
# This does not actually allocate
function get_free_slot(mem::CodeMem, allow_inplace = [])
    for i = 1:size(mem.slots, 1)
        if mem.slots[i] == :Free || mem.slots[i] in allow_inplace
            return (i, mem.slot_names[i])
        end
    end
    return error("No free memory slot found!")
end

# Allocates a slot that has been found using
# get_free_slot
function alloc_slot!(mem::CodeMem, i, node::Symbol)
    return mem.slots[i] = node
end

function free!(mem::CodeMem, i)
    return mem.slots[i] = :Free
end
