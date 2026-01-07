"""
    create_return_time_plot(vals, names, limits, title, filename)

Generates a custom "return time measure" plot, which visualizes samples 
against a defined range (min, expected, max) on a vertical line.

# Arguments
- `vals::Vector{<:Number}`: The y-coordinates (sample values) to plot.
- `names::Vector{<:AbstractString}`: The text labels for each sample in `vals`.
- `limits::Vector{<:Number}`: A vector containing [min_val, expected_val, max_val] in that order.
- `title::AbstractString`: The main title for the plot.
- `filename::AbstractString`: The path (e.g., "results/my_plot.png") to save the figure to.

# Returns
- A `Plot` object.
"""
# function plot_measure_nodes(vals::Vector{<:Number},
#                                  names::Vector{<:AbstractString},
#                                  limits::Vector{<:Number},
#                                  title::AbstractString,
#                                  filename::AbstractString)
    
#     # --- 1. Input Validation and Extraction ---
#     if length(limits) != 3
#         error("The 'limits' vector must contain exactly three values: [min_val, expected_val, max_val].")
#     end
#     min_val, expected_val, max_val = limits

#     # Constants used in the original plot styling
#     tick_width = 0.03
    
#     # --- 2. Base Plot and Vertical Line ---
#     # We store the initial plot in `p` so subsequent calls use `!` (mutating functions)
#     p = plot([0, 0], [min_val, max_val];
#              color=:black, linewidth=2,
#              # Using the original hardcoded limits for the axis padding
#              ylim=(min_val - 0.1, max_val + 0.2), 
#              xlim=(-0.6, 0.6),
#              legend=false, xlabel="", ylabel="", framestyle=:none,
#             #  title=title,
#              size=(500, 800))

#     # --- 3. Expected Value Marker ---
#     scatter!([0], [expected_val];
#              color=:white, markerstrokecolor=:black,
#              marker=:circle, markersize=11, markerstrokewidth=2, 
#              label="Expected", 
#              parent=p) # Explicitly drawing onto the plot p

#     # --- 4. Sample Values and Names ---
#     # X-coordinates for sample points (shifted right from the center line)
#     xvals = fill(0.08, length(vals))
#     scatter!(xvals, vals; color=:red, marker=:diamond, label="Samples", parent=p)

#     # # Horizontal annotations for the sample names (to the right)
#     # for (x, y, name) in zip(xvals, vals, names)
#     #     # Shift annotation slightly further right
#     #     annotate!(x + 0.04, y, text(name, 8, :left), parent=p) 
#     # end

#     # --- 5. Range Boundary Ticks ---
#     # Small horizontal lines at min and max
#     plot!([-tick_width, tick_width + 0.05], [min_val, min_val]; color=:black, linewidth=2, parent=p)
#     plot!([-tick_width, tick_width + 0.05], [max_val, max_val]; color=:black, linewidth=2, parent=p)

#     # --- 6. Numeric and Semantic Labels (Left Side) ---
    
#     # Numeric labels (rounded to 2 digits)
#     annotate!(-0.05, min_val, text(string(round(min_val, digits=2)), 8, :right), parent=p)
#     annotate!(-0.05, expected_val, text(string(round(expected_val, digits=2)), 8, :right), parent=p)
#     annotate!(-0.05, max_val, text(string(round(max_val, digits=2)), 8, :right), parent=p)

#     # Semantic labels 
#     annotate!(-0.2, min_val, text("min", 8, :right), parent=p)
#     annotate!(-0.2, expected_val, text("expected", 8, :right), parent=p)
#     annotate!(-0.2, max_val, text("max", 8, :right), parent=p)

#     # --- 7. Saving the Figure ---
#     savefig(p, filename)
    
#     println("Plot saved to: $filename")

#     return p
# end



# takes normalized values and makes the plot
function plot_norm_SensInf(taxa_names, sens, infl, filename)
    plt = groupedbar(
        # taxa_names,
        [sens infl],
        label = ["sensitivity" "influence"],
        xticks=(1:length(taxa_names), taxa_names),  # taxa names on x-axis
        xrotation=45,
        size=(1000, 600),
        # yticks=([-1.0, 0.0, 1.0], ["not at all", "average", "the most possible"]),
        ylabel = "normalized over-/under-average",
        ylimits = [-1, 1],
        #  title = "Sensitivity and Influence",
        bar_width = 0.27,
        grouped = true
        # xmirror=true,   # mirror x-axis to the top
        # orientation = :horizontal   # make bars horizontal
    )
    hline!([0.0], color=:black, lw=2, label=false)

    # --- 7. Saving the Figure ---
    savefig(plt, filename)
    
    println("Plot saved to: $filename")
    return plt
end


function make_barplot(values::AbstractVector, xticks::Vector{String}, ylabel::String, filename::String, limits::Vector{Float64}= Float64[])
    if length(values) != length(xticks)
        error("Length of values and xticks must match.")
    end

    # plt = bar(
    #     values,
    #     xticks = (1:length(values), xticks),
    #     xrotation=45,
    #     ylabel = ylabel,
    #     legend = false,
    #     color = :gray,
    #     bar_width = 0.3,
    # )
    n = length(values)
    colors = Vector{Any}(undef, n)
    markers = fill(:circle, n)
    markercolors = Vector{Any}(undef, n)
    markerstrokecolors = fill(:black, n)

    # #Arctic
    # colors[1]        = :black #detritus
    # colors[2:3]      .= :gray50 #bacteria
    # colors[4:10]     .= :green
    # colors[11:37]    .= :blue #mesozooplankton
    # colors[38:39]    .= :gray80 #fish
    # colors[40:41]    .= :brown #seabirds
    # colors[42:45]    .= :black #mammals
    # markercolors[1:41] .= colors[1:41]      # normal fill
    # markercolors[42:45] .= :white           # white fill
    # markerstrokecolors[42:45] .= :black

    #GreenlandSea
    colors[1]        = :black
    colors[2:3]      .= :gray50
    colors[4:10]     .= :green
    colors[11:38]    .= :blue
    colors[39:45]    .= :gray80
    colors[46:48]    .= :brown
    colors[49:56]    .= :gray20
    markercolors[1:48] .= colors[1:48]      # normal fill
    markercolors[49:56] .= :white           # white fill
    markerstrokecolors[49:56] .= :black

    plt = bar(
        1:length(values),        # x positions
        values,                  # y values
        xticks = (1:length(values), xticks),
        xtickfontsize = 11,
        xrotation = 45,
        ylabel = ylabel,
        yscale =:log10,
        legend = false,
        color = colors,
        # marker = (:circle, 4),   # circle markers, size 8
        # marker = (markers, 4),
        # markercolor = markercolors,
        # markerstrokecolor = markerstrokecolors,
        size = (n*20+280, 600),   # ← width, height in pixels
    )
    plot!(plt, xtickfontsize = 11)
#     plt = bar(
#     values,                          # bar lengths (x)
#     1:length(values),                # y positions
#     yticks = (1:length(values), ticks),
#     ylabel = "",                     # labels now on y-axis
#     xlabel = label,
#     xscale = :log10,
#     legend = false,
#     orientation = :horizontal,
#     color = colors,
# )
    # Add horizontal limit lines if provided
    if !isempty(limits)
        if length(limits) != 3
            error("limits must be a vector of length 3: [min_val, exp_val, max_val]")
        end
        
        min_val, exp_val, max_val = limits

        hline!(plt, [min_val], linestyle=:solid, linewidth=2, color=:black,label="min")
        hline!(plt, [exp_val], linestyle=:dash, linewidth=2, color=:black,label="expected")
        hline!(plt, [max_val], linestyle=:solid, linewidth=2, color=:black, label="max")

        # Extract current yticks safely
        ticks_tuple = yticks(plt)[1]
        current_yticks = ticks_tuple[1]  # numeric positions
        # println(current_yticks)
        current_labels = length(ticks_tuple) > 1 ? ticks_tuple[2] : string.(ticks_tuple[1])

        # Combine current ticks with limit positions
        combined_ticks = sort(union(current_yticks, [min_val, exp_val, max_val]))
        
        # Build corresponding labels: default to numeric label, replace limit positions with strings
        combined_labels = string.(combined_ticks)  # start with numeric labels
        for (i, t) in enumerate(combined_ticks)
            if t == min_val
                combined_labels[i] = "min"
            elseif t == exp_val
                combined_labels[i] = "expected"
            elseif t == max_val
                combined_labels[i] = "max"
            end
        end

        # Apply new yticks
        yticks!(plt, combined_ticks, combined_labels)
    end

    # # Determine lower and upper bounds
    # data_min = minimum(values)
    # data_max = maximum(values)

    # if !isempty(limits)
    #     data_min = min(data_min, minimum(limits))
    #     data_max = max(data_max, maximum(limits))
    # end

    # # Padding factors (log-scale friendly)
    # lower = data_min / 1.5
    # upper = data_max * 1.5

    # # Apply y-limits
    # ylims!(plt, lower, upper)

    # # Make bars start from the bottom of the plot
    # bar!(plt; baseline = lower)

    savefig(plt, filename)
    
    println("Plot saved to: $filename")
    return plt
end

function plot_reactivity(values::AbstractVector, xticks::Vector{String}, ylabel::String, filename::String, limits::Vector{Float64}, resilience)
    if length(values) != length(xticks)
        error("Length of values and xticks must match.")
    end

    # plt = bar(
    #     values,
    #     xticks = (1:length(values), xticks),
    #     xrotation=45,
    #     ylabel = ylabel,
    #     legend = false,
    #     color = :gray,
    #     bar_width = 0.3,
    # )
    plt = scatter(
        1:length(values),        # x positions
        values,                  # y values
        xticks = (1:length(values), xticks),
        xrotation = 45,
        ylabel = ylabel,
        # yscale =:log10,
        legend = false,
        color = :gray,
        marker = (:circle, 4),   # circle markers, size 8
    )
    # Add horizontal limit lines if provided
    
        
        min_val, exp_val, max_val = limits

        hline!(plt, [min_val], linestyle=:solid, linewidth=2, color=:black,label="min")
        hline!(plt, [exp_val], linestyle=:dash, linewidth=2, color=:black,label="expected")
        hline!(plt, [max_val], linestyle=:solid, linewidth=2, color=:black, label="reactivity")
        hline!(plt, [resilience], linestyle=:dot, linewidth=2, color=:black, label="-resilience")

        # Extract current yticks safely
        ticks_tuple = yticks(plt)[1]
        current_yticks = ticks_tuple[1]  # numeric positions
        # println(current_yticks)
        current_labels = length(ticks_tuple) > 1 ? ticks_tuple[2] : string.(ticks_tuple[1])

        # Combine current ticks with limit positions
        combined_ticks = sort(union(current_yticks, [min_val, exp_val, max_val, resilience]))
        
        # Build corresponding labels: default to numeric label, replace limit positions with strings
        combined_labels = string.(combined_ticks)  # start with numeric labels
        for (i, t) in enumerate(combined_ticks)
            if t == min_val
                combined_labels[i] = "min"
            elseif t == exp_val
                combined_labels[i] = "expected"
            elseif t == max_val
                combined_labels[i] = "reactivity"
            elseif t == resilience
                combined_labels[i] = "-resilience"
            end
        end

        # Apply new yticks
        yticks!(plt, combined_ticks, combined_labels)
    

    savefig(plt, filename)
    
    println("Plot saved to: $filename")
    return plt
end

function compare_comm_plot(ymin::AbstractVector, ymax::AbstractVector, xticks::Vector{String}, ylabel::String, filename::String, med::Vector{Float64}= Float64[])
    # choose cap width in x-data units (10% of minimum x spacing)
    # dx = minimum(diff(sort(x)))
    cap = 0.03    # tweak 0.08..0.2 as needed

    plt = plot(xlim=(minimum(x)-1, maximum(x)+1), ylim=(minimum(ymin), maximum(ymax)))
    for (xi, ylo, yhi) in zip(x, ymin, ymax)
        # vertical line
        plot!(plt, [xi, xi], [ylo, yhi], color=:black, linewidth=2, label=false)
        # bottom cap
        plot!(plt, [xi - cap, xi + cap], [ylo, ylo], color=:black, linewidth=2, label=false)
        # top cap
        plot!(plt, [xi - cap, xi + cap], [yhi, yhi], color=:black, linewidth=2, label=false)
    end

    # optional: add marker for median/mean
    if !isempty(med)
        scatter!(plt, x, med, color=:black, marker=:circle, markersize = 6, label=false)
    end
    # rotated labels

    plot!(
        plt,
        xticks = (x, xticks),
        xrotation = 45,     # optional: rotate labels if long
        yscale = :log10,       # << log10 scale
        ylabel = ylabel
    )


    savefig(plt, filename)
    
    println("Plot saved to: $filename")
    return plt
end
