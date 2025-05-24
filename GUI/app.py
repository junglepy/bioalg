from flask import Flask, render_template, request

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def index():
    result = None
    if request.method == "POST":
        try:
            number = int(request.form["number"])
            result = number + 1
        except ValueError:
            result = "Ошибка: введите целое число"
    return render_template("index.html", result=result)

if __name__ == "__main__":
    app.run(port=7700)
