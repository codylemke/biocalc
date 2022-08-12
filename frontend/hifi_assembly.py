"""Reaction mix calculator for PCR

Description...

Usage Example:
    ...Example
"""
# Third Party Packages
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import pandas as pd
# Global Constants
app = Dash(__name__)


# LAYOUT ----------------------------------------------------------------------
app.layout = html.Div([
    html.H2('HiFi Assembly'),
    dash_table.DataTable(
        id='parameters',
        columns=(
            [
                {'id': 'Key'}
            ]
        ),
        data=[
            {'Key': 'Value'}
        ],
        editable=True,
        style_data_conditional=[
            {
                'if': {'column_editable': True},
                'backgroundColor': 'rgb(225, 225, 50)',
            }
        ],
    ),
    html.Br(),
    html.H2('Entries'),
    dash_table.DataTable(
        id='entries',
        columns=[
            {'id': 'Key'},
        ],
        data=[
            {'id': 'Key'}
        ],
    ),
])


# CALLBACKS -------------------------------------------------------------------
@app.callback(
    Output('parameters', 'data'),
    Input('table-editing-simple', 'data'),
    Input('table-editing-simple', 'columns'))
def update_table(rows, columns):
    """Placeholder"""
    return


# SCRIPT ----------------------------------------------------------------------
if __name__ == '__main__':
    app.run(host='127.0.0.1', debug=True)
