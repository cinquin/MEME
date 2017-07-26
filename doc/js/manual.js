function toggle_man_display() {
  var index;
  var ids = ["man_web_button", "man_cmd_button", "man_usage", "man_cmd"];
  for (index = 0; index < ids.length; ++index) {
    var element = document.getElementById(ids[index]);
      element.style.display=(element.style.display=='none') ? element.style.display='inline' : element.style.display='none';
  }
}

(function () {
  var parts = window.location.search.substr(1).split("&");
  var $_GET = {};
  for (var i = 0; i < parts.length; i++) {
      var temp = parts[i].split("=");
      $_GET[decodeURIComponent(temp[0])] = decodeURIComponent(temp[1]);
  }
  var man_type = $_GET['man_type'];
  //alert(man_type);
  if (man_type == 'cmd') {
    man_cmd_button.style.display='none'
    man_web_button.style.display='inline'
    man_cmd.style.display='inline'; 
    man_usage.style.display='inline'
  }
  if (man_type == 'web') {
    man_cmd_button.style.display='inline'
    man_web_button.style.display='none'
    man_cmd.style.display='none'
    man_usage.style.display='none'
  }
})();
