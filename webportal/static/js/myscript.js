
// // SESSION CREATE JAVASCRIPT
$( document ).ready(function() {
  var user_url = window.location.pathname
  if (user_url.includes('session_create')){ // executes if session_create page loaded
    console.log('session_create page loaded')
    $('#id_select_genome').hide()
    $('#id_organism').hide()
    $('#id_salmon').hide()
    $('#id_fasta_dna_file').hide()
    $('#id_fasta_cdna_file').hide()
    $('#id_gtf_file').hide()
    $('#id_genome_index').on('click', function(){
      console.log('option selected')
      if ($('#id_genome_index option:selected').val() == 'pre_index'){ // executes if pre_index selected from dropdown
        console.log('great success with pre index')
        $('#id_organism').hide()
        $('#id_salmon').hide()
        $('#id_fasta_dna_file').hide()
        $('#id_fasta_cdna_file').hide()
        $('#id_gtf_file').hide()
        $('#id_select_genome').show()
      }
      if ($('#id_genome_index option:selected').val() == 'user_provided'){ // executes if pre_index selected from dropdown
        console.log('great success with user provided')
        $('#id_select_genome').hide()
        $('#id_organism').show()
        $('#id_salmon').show()
        $('#id_salmon').on('click', function(){
          if ($('#id_salmon option:selected').val() == 2){
            console.log('salmon selected')
            $('#id_fasta_dna_file').show()
            $('#id_fasta_cdna_file').show()
            $('#id_gtf_file').show()
          }
          if ($('#id_salmon option:selected').val() == 3){
            console.log('salmon NOT selected')
            $('#id_fasta_cdna_file').hide()
            $('#id_fasta_dna_file').show()
            $('#id_gtf_file').show()
          }
        })
      }
    })
  }
})




if ($('#div_id_libtype option:selected').val() == 'user_provided'){ // executes if user_provided selected from dropdown
  $('#div_id_read_1').show()
  $('#div_id_accession').show()
  $('#div_id_read_2').hide()
}


// SAMPLES CREATE JAVASCRIPT
$( document ).ready(function() {
  var user_url = window.location.pathname
  if (user_url.includes('samples_create')){ // executes if samples_create page loaded
    $('#div_id_read_1').hide()
    $('#div_id_read_2').hide()
    $('#div_id_accession').hide()
    $('#validate_library').hide()
    $('#div_id_libtype').on('click', function(){
      if ($('#div_id_libtype option:selected').val() == 'PE'){ // executes if PE selected from dropdown
        $('#div_id_read_1').show()
        $('#div_id_read_2').show()
        $('#div_id_accession').show()
      }
      if ($('#div_id_libtype option:selected').val() == 'SG'){ // executes if SG selected from dropdown
        $('#div_id_read_1').show()
        $('#div_id_accession').show()
        $('#div_id_read_2').hide()
      }
    })
  }
});



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
