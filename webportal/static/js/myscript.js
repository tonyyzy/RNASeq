

// FUNCTIONS
var i;
function unhide(divID) {
  console.log('unhiding')
  for (i = 0; i < divID.length; i++){
    div = divID[i]
    console.log(div)
    var item = document.getElementById(div);
    if (item) {
        item.className = 'unhidden' ;
    }
  }
};
function hide(divID) {
  console.log('hiding')
  for (i = 0; i < divID.length; i++){
    div = divID[i]
    console.log(div)
    var item = document.getElementById(div);
    if (item) {
        item.className = 'hidden';
    }
  }
}


// SAMPLES CREATE JAVASCRIPT
// $( document ).ready(function() {
//   var user_url = window.location.pathname
//   if (user_url.includes('samples_create')){ // executes if samples_create page loaded
//     $('#div_id_read_1').hide()
//     $('#div_id_read_2').hide()
//     $('#div_id_accession').hide()
//     $('#validate_library').hide()
//     $('#div_id_libtype').on('click', function(){
//       if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
//         $('#div_id_read_1').show()
//         $('#div_id_read_2').show()
//         $('#div_id_accession').show()
//       }
//       if ($('#div_id_libtype option:selected').val() == 'SG'){ // executes if SG selected from dropdown
//         $('#div_id_read_1').show()
//         $('#div_id_accession').show()
//         $('#div_id_read_2').hide()
//       }
//     })
//   }
// });



// SAMPLES UPDATE JAVASCRIPT
$( document ).ready(function() {
  var user_url = window.location.pathname // stores current url as user_url
  if (user_url.includes('samples_update')){ // executes if samples_update page loaded
    $('#validate_library').hide()
    if ($('#div_id_libtype option:selected').val() == 'SG'){
    	$('#div_id_read_2').toggle()
    }
    $('#div_id_libtype').on('click', function(){
      if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
        $('#div_id_read_2').show()
        $('#read_2-clear_id').prop("checked", false);
      }
      if ($('#div_id_libtype option:selected').val() == 'SG'){ // executes if SG selected from dropdown
        $('#div_id_read_2').hide()
        $('#read_2-clear_id').prop("checked", true);
      }
    })
    }
});
