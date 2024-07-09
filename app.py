from flask import Flask, render_template


app = Flask(__name__)

JOBS = [
    { 'title': 'Software Engineer', 'company': 'Google', 'location': 'Mountain View, CA', 'date_posted': '2021-01-01', 'description': 'We are looking for a software engineer to join our team.'},
    { 'title': 'Data Scientist', 'company': 'Facebook', 'date_posted': '2021-01-02', 'description': 'We are looking for a data scientist to join our team.'},
    { 'title': 'Product Manager', 'company': 'Apple', 'location': 'Cupertino, CA', 'date_posted': '2021-01-03', 'description': 'We are looking for a product manager to join our team.'},
]

@app.route('/')

def hello():
    return  render_template('home.html', jobs=JOBS)

@app.route('/home')
def home():
    return 'Welcome to the Home Page!'



if __name__ == '__main__':
    app.run(debug=True)