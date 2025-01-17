from os.path import dirname, join

import pandas as pd
import numpy as np
import colorcet as cc
import random

from bokeh.plotting import figure
from bokeh.layouts import layout, column
from bokeh.models import ColumnDataSource, Div
from bokeh.models.widgets import Slider, Select, TextInput
from bokeh.io import curdoc

#movies["color"] = np.where(movies["Oscars"] > 0, "orange", "grey")
#movies["alpha"] = np.where(movies["Oscars"] > 0, 0.9, 0.25)
#movies.fillna(0, inplace=True)  # just replace missing values with zero
#movies["revenue"] = movies.BoxOffice.apply(lambda x: '{:,d}'.format(int(x)))

tsne = pd.read_csv(join(dirname(__file__), "../results/tsne_diag.csv"))
#movies.loc[movies.imdbID.isin(razzies), "color"] = "purple"
#movies.loc[movies.imdbID.isin(razzies), "alpha"] = 0.9

axis_map = {
    "t-SNE 1": "tsne1",
    "t-SNE 2": "tsne2",
#    "Number of Reviews": "Reviews",
#    "Box Office (dollars)": "BoxOffice",
#    "Length (minutes)": "Runtime",
#    "Year": "Year",
}

random.seed(30)
meth_grp_uniq = list(set(tsne['grp']))
colors = random.sample(cc.glasbey_dark, len(meth_grp_uniq))

color_assign = []
alpha_assign = []
for group in tsne['grp']:
    if group == "DIAGNOSTIC":
        color_assign.append("black")
        alpha_assign.append(1.0)
    else:
        color_assign.append(colors[meth_grp_uniq.index(group)])
        alpha_assign.append(0.5)

tsne['color'] = color_assign
tsne['alpha'] = alpha_assign

desc = Div(text=open(join(dirname(__file__), "description.html")).read(), sizing_mode="stretch_width")

# Create Input controls
reviews = Slider(title="Minimum number of reviews", value=80, start=10, end=300, step=10)
min_year = Slider(title="Year released", start=1940, end=2014, value=1970, step=1)
max_year = Slider(title="End Year released", start=1940, end=2014, value=2014, step=1)
oscars = Slider(title="Minimum number of Oscar wins", start=0, end=4, value=0, step=1)
boxoffice = Slider(title="Dollars at Box Office (millions)", start=0, end=800, value=0, step=1)
genre = Select(title="Genre", value="DIAGNOSTIC",
               options=meth_grp_uniq + ["All"])
director = TextInput(title="Director name contains")
cast = TextInput(title="Cast names contains")
x_axis = Select(title="X Axis", options=sorted(axis_map.keys()), value="t-SNE 1")
y_axis = Select(title="Y Axis", options=sorted(axis_map.keys()), value="t-SNE 2")

# Create Column Data Source that will be used by the plot
source = ColumnDataSource(data=dict(x=[], y=[], color=[], grp=[], alpha=[]))

#TOOLTIPS=[
#    ("Title", "@title"),
#    ("Year", "@year"),
#    ("$", "@revenue")
#]

#p = figure(plot_height=600, plot_width=700, title="", toolbar_location=None, tooltips=TOOLTIPS, sizing_mode="scale_both")
p = figure(plot_height=600, plot_width=700, title="", toolbar_location=None, sizing_mode="fixed")
p.circle(x="x", y="y", source=source, size=7, color="color", line_color=None, fill_alpha="alpha")


def select_movies():
    genre_val = genre.value
#    director_val = director.value.strip()
#    cast_val = cast.value.strip()
    selected = tsne#[
#        (movies.Reviews >= reviews.value) &
#        (movies.BoxOffice >= (boxoffice.value * 1e6)) &
#        (movies.Year >= min_year.value) &
#        (movies.Year <= max_year.value) &
#        (movies.Oscars >= oscars.value)
#    ]
    if (genre_val != "All"):
        selected = selected[selected.grp.str.contains(genre_val)==True]
#    if (director_val != ""):
#        selected = selected[selected.Director.str.contains(director_val)==True]
#    if (cast_val != ""):
#        selected = selected[selected.Cast.str.contains(cast_val)==True]
    return selected


def update():
    df = select_movies()
    x_name = axis_map[x_axis.value]
    y_name = axis_map[y_axis.value]

    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value
    p.title.text = "%d movies selected" % len(df)
    source.data = dict(
        x=df[x_name],
        y=df[y_name],
        color=df["color"],
#        title=df["Title"],
#        year=df["Year"],
#        revenue=df["revenue"],
        alpha=df["alpha"],
    )

controls = [reviews, boxoffice, genre, min_year, max_year, oscars, director, cast, x_axis, y_axis]
for control in controls:
    control.on_change('value', lambda attr, old, new: update())

inputs = column(*controls, width=320, height=1000)
inputs.sizing_mode = "fixed"
l = layout([
    [desc],
    [inputs, p],
], sizing_mode="scale_both")

update()  # initial load of the data

curdoc().add_root(l)
curdoc().title = "Movies"
