"""Generate Golden Gate Express Fragments

Descriptions

Usage Example:
    example
"""
# Third Party Packages
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output, State
# Local Modules
from models import DNA
# Global Constants
app = Dash(__name__)


# LAYOUT ----------------------------------------------------------------------
app.layout = html.Div([
    html.H1('Golden Gate Express Fragment Generator'),
    html.P('Number of Sequences'),
    dcc.Input(
        id='number-of-sequences',
        type='number',
        placeholder='number_of_sequences'
    ),
    html.Br(),
    html.H2('Inputs'),
    html.Button('Submit Sequences', id='submit-sequences', n_clicks=0),
    dash_table.DataTable(
        id='input-table',
        columns=[
            {'name': 'Name', 'id': 'name', 'type': 'text', 'clearable': True},
            {'name': 'Organism', 'id': 'organism', 'type': 'text', 'clearable': True},
            {'name': 'Module', 'id': 'module', 'type': 'numeric', 'clearable': True},
            {'name': 'Sequence', 'id': 'sequence', 'type': 'text', 'clearable': True},
        ],
        data=[
            {'name': '', 'organism': '', 'module': '', 'sequence': ''}
            for num in range(10)
        ],
        editable=True,
    ),
    html.H2('Output Sequences'),
    dash_table.DataTable(
        id='output-table',
        columns=[
            {'name': 'Name', 'id': 'name', 'type': 'text'},
            {'name': 'Organism', 'id': 'organism', 'type': 'text'},
            {'name': 'Module', 'id': 'module', 'type': 'numeric'},
            {'name': 'GGE Sequence', 'id': 'gge_sequence', 'type': 'text'},
        ],
        column_selectable='single',
        data=[]
    ),
])


# CALLBACKS -------------------------------------------------------------------
@app.callback(
    Output('input-table', 'data'),
    Input('number-of-sequences', 'value'),)
def update_table(value):
    """Placeholder"""
    return [
        {'name': '', 'module': '', 'organism': '', 'sequence': ''}
        for num in range(value)
    ]

@app.callback(
    Output('output-table', 'data'),
    Input('submit-sequences', 'n_clicks'),
    State('input-table', 'data'))
def process_request(n_clicks, data):
    """Runs all of the sequences through append adapters"""
    new_data = list()
    for row in data:
        dna = DNA(name=row['name'], sequence=row['sequence'])
        fragment = dna.generate_gge_fragment(organism=row['organism'], module=row['module'])
        new_data.append({
            'name': row['name'],
            'organism': row['organism'],
            'module': row['module'],
            'gge_sequence': fragment.sequence
        })
    return new_data


# SCRIPT ----------------------------------------------------------------------
if __name__ == '__main__':
    app.run(host='127.0.0.1', debug=True)
