
from importlib.resources import files

import dash
from dash import dcc, html, dash_table, Input, Output, State, ctx, callback
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from cmd_fusion.data_objects.gene_set import GeneSet
from cmd_fusion.dashboard.arguments import parse_dashboard_args
from cmd_fusion.dashboard.general import head_layout
from cmd_fusion.dashboard.data import (_read_default_cols, _load_json_table,
                                       _retrieve_cases, _read_display_names, _read_case)


args = parse_dashboard_args()
case_dir = args["sample_dir"]

print("Reading case files...")
case_list = _retrieve_cases(case_dir)

print("Reading gene set...")
if ".pkl" in args["gtf"]:
    geneset = GeneSet.read_pickle(args["gtf"])
else:
    geneset = GeneSet(args["gtf"])

print("Reading column data...")
display_names = _read_display_names()
default_cols = _read_default_cols()
tab_dict = {"overview_tab": default_cols,
            "thermo_tab": ["General", "Thermo"],
            "arriba_tab": ["General", "Arriba"],
            "star_tab": ["General", "STAR"]}


print("Building app...")
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, "assets/tabs.css"])
app.title = "CMD Fusion Dashboard"
app.layout = html.Div([
    head_layout,
    html.Div(html.H1("Case Explorer")),
    dcc.Download(id="download-content"),

    dbc.Row(justify="start", children=[
            dbc.Col(sm=6, md=4, lg=3, xl=2, children=[
                dcc.Dropdown(case_list, id="case-dropdown", placeholder="Select a case")
                    ]),
            dbc.Col(children=[
                html.H3("Case ID:")
                ]),
            dbc.Col(children=[
                dcc.Loading(children=[
                    html.Div([html.H3("None")], id="case-id",
                             loading_state=dict(component_name="display-table"))
                    ])
                ]),
            dbc.Col(children=[
                html.Button("Download View - TSV", id="download-view", n_clicks=0),
                html.Button("Download All - TSV", id="download-all", n_clicks=0)
                ])
            ]),

    html.Br(),

    dcc.Tabs(id="tool-tabs",
             parent_className="custom-tabs",
             className="custom-tabs-container",
             mobile_breakpoint=0,
             value=None,
             children=[
                 dcc.Tab(id=f"{seq}_tab",
                         className="custom-tab",
                         selected_className="custom-tab--selected",
                         label=seq.upper(),
                         value=f"{seq}_tab",
                         disabled=False,
                         children=[])
                 for seq in ["overview", "thermo", "arriba", "star"]]),

    html.Div(id="display-table"),
    ])


@callback(
    Output("current-table", "data"),
    Output("tool-tabs", "value"),
    Output("case-id", "children"),
    Input("case-dropdown", "value")
    )
def load_case(case):
    if case is None:
        raise PreventUpdate
    case_id = html.H3(case)
    case = _read_case(case, case_dir, geneset)
    return case, "overview_tab", case_id


def prep_dash_table(table, tool):
    columns = tab_dict[tool]
    displayed = table[[col for col in columns if col in table.columns]]
    displayed = displayed[[col for col in display_names.keys()
                           if col in displayed.columns]]
    cols = [{"name": display_names[col], "id": ",".join(col)}
            for col in displayed.columns]
    for col in cols:
        if col["id"] == "Thermo,Annotation,COSMIC":
            col["presentation"] = "markdown"

    # if tool != "overview_tab":
    #     matches = displayed[("General", "FoundBy")].str.contains(columns[1])
    #     displayed = displayed.loc[matches]

    displayed.columns = [",".join(col) for col in displayed.columns]
    displayed = displayed.to_dict("records")
    prepped = {"data": displayed, "columns": cols}
    return prepped


@callback(
    Output("display-table", "children"),
    Input("tool-tabs", "value"),
    State("current-table", "data")
    )
def choose_tool(tool, table_json):
    table = _load_json_table(table_json)
    prepped = prep_dash_table(table, tool)
    table = dash_table.DataTable(
        **prepped,
        style_table={"overflowX": "auto"},
        markdown_options={"html": True},
        merge_duplicate_headers=True,
        page_action="native", page_size=50,
        filter_action="native",
        filter_options={"case": "insensitive"},
        sort_action="native",
        )
    return table


@callback(
    Output("download-content", "data"),
    Input("download-view", "n_clicks"),
    Input("download-all", "n_clicks"),
    State("current-table", "data"),
    State("tool-tabs", "value"),
    State("case-id", "children"),
    prevent_initial_call=True
    )
def download_view(clicks1, clicks2, table_json, tool, case_id):
    table = _load_json_table(table_json)
    tool_str = ""
    if ctx.triggered_id == "download-view":
        cols = [col for col in tab_dict[tool] if col in table.columns]
        table = table[cols]
        tool_str = "." + tool.split('_')[0]
    case_id = case_id["props"]["children"]
    return dcc.send_data_frame(table.to_csv,
                               f"{case_id}.fusion_dash{tool_str}.tsv",
                               sep="\t", index=True)

print("Running app...")
app.run_server()
