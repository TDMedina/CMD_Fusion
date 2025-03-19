
import dash_bootstrap_components as dbc
from dash import dcc, html

memory = html.Div(children=[dcc.Store(id="current-table", data={})])

nav_title = dcc.Link(href="home",
                     style={"textDecoration": "none"},
                     children=[
                         dbc.Row(align="center",
                                 className="g-0",
                                 children=[
                                     dbc.Col(children=[
                                         # html.Img(src="assets/sjcmd.png", height="50px")
                                         ]),
                                     dbc.Col(children=[
                                         dbc.NavbarBrand("CMD Fusion Dashboard",
                                                         className="ms-2")
                                         ])
                                     ])
                         ])

navbar = dbc.Navbar(
    children=[dbc.Container(
        children=[nav_title]
        )],
    color="#1B2F52",
    dark=True,
    className="mb-5",
    )

head_layout = html.Div([memory,
                        dcc.Location(id="url", refresh=False),
                        navbar])
