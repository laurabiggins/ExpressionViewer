$(document).on("shiny:inputchanged", function(event) {
  if (event.inputType === "shiny.number" && (event.name === "WidthPixels" || event.name === "HeightPixels")) {
    if (event.value === null || event.value < 100 || event.value > 10000) {
      $("#" + event.name).parents(".form-group").addClass("has-error");
    } else {
      $("#" + event.name).parents(".form-group").removeClass("has-error");
    }
  }

  if (event.inputType === "shiny.number" && event.name === "log2FCFilt") {
    if (event.value === null || event.value < 0 || event.value > 20) {
      $("#" + event.name).parents(".form-group").addClass("has-error");
    } else {
      $("#" + event.name).parents(".form-group").removeClass("has-error");
    }
  }

  if (event.inputType === "shiny.number" && event.name === "AdjpFilt") {
    if (event.value === null || event.value < 0 || event.value > 1) {
      $("#" + event.name).parents(".form-group").addClass("has-error");
    } else {
      $("#" + event.name).parents(".form-group").removeClass("has-error");
    }
  }
});