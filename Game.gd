extends Node2D

var bigG=  3000
var display_step = 0.1

func _ready():
	$CanvasLayer/labelTimeScale.text = str(Engine.time_scale)
	
	#$KPlanet.eccentricty=0.99
	#$KPlanet.guessE(0.19006693557797, 0.19006693557797)
	
	#$KPlanet.eccentricty=0.99
	#$KPlanet.guessE(0.19189450226622, 0.19189450226622)
	#$KPlanet.guessE(0.12427453480098, 0.12427453480098)
	
	#$KPlanet.eccentricty=0.98
	#$KPlanet.guessE(0.14986046843648, 0.14986046843648)
	
var totalDelta	= 0.0	
var gameTick = 0
func _process(delta):
	totalDelta += delta
	if gameTick < int(totalDelta):
		gameTick =  int(totalDelta)
	$CanvasLayer/labelTotalDelta.text = str(gameTick)
	#$Camera2D.position = get_node_or_null('Star').global_position

var bypass_slider_signals = false
func _on_HSliderSemiMajorAxis_value_changed(value):
	if (bypass_slider_signals):
		return	
	$KPlanet.semiMajorAxis = value
	$KPlanet.solveFromTrue()
	$KPlanet.queue_redraw()
func _on_HSliderTrue_value_changed(value):
	if (bypass_slider_signals):
		return
	$KPlanet.trueAnomaly = value
	$KPlanet.solveFromTrue()
	$KPlanet.queue_redraw()
func _on_HSliderArgPeri_value_changed(value):
	if (bypass_slider_signals):
		return	
	$KPlanet.argPeri = value
	$KPlanet.solveFromTrue()	
	$KPlanet.queue_redraw()
func _on_HSliderEccent_value_changed(value):
	if (bypass_slider_signals):
		return	
	$KPlanet.eccentricty=value
	$KPlanet.solveFromTrue()	
	$KPlanet.queue_redraw()	
func _on_HSliderEccAnom_value_changed(value):
	if (bypass_slider_signals):
		return		
	$KPlanet.eccentricAnomaly = value
	$KPlanet.solveFromEccentricAnom()	
	$KPlanet.queue_redraw()
func _on_HSliderMeanAnom_value_changed(value):
	if (bypass_slider_signals):
		return		
	$KPlanet.meanAnomaly = value
	$KPlanet.solveFromMeanAnom()	
	$KPlanet.queue_redraw()
func _on_HSliderBigG_value_changed(value):
	if (bypass_slider_signals):
		return		
	$KPlanet.setG(value)	
	$KPlanet.solveFromTrue()	
	$KPlanet.queue_redraw()
func _on_HSlider_value_changed(value):
	$CanvasLayer/labelTimeScale.text = str(value)
	Engine.time_scale = value

func _unhandled_input(event):
	if event is InputEventMouseButton:
		if event.button_index == MOUSE_BUTTON_WHEEL_UP and event.pressed:
			zoomIn()
		elif event.button_index == MOUSE_BUTTON_WHEEL_DOWN and event.pressed:
			zoomOut()
		   
func zoomIn():
#	var _zoom_pos 6= get_global_mouse_position()
	if snapped($Camera2D.zoom.x, .1) < .6:
		$Camera2D.zoom = $Camera2D.zoom + (Vector2.ONE/10)

func zoomOut():
#	var _zoom_pos = get_global_mouse_position()
	if snapped($Camera2D.zoom.x, .1) > .1:
		$Camera2D.zoom = $Camera2D.zoom - (Vector2.ONE/10)

var logMsgs = Array()
var logMaxMsgs = 10
func logMsg(msg,clear_first=false):
	if clear_first:
		logMsgs.clear()
	logMsgs.append(msg)
	if logMsgs.size() > logMaxMsgs:
		logMsgs.remove_at(0)
	$CanvasLayer/LogPanel/labelLog.text = "\n".join(PackedStringArray(logMsgs))

func _on_Game_draw():
	draw_arc($Star.global_position, $Star.max_distance_allowed, 0, TAU, 1000, Color.WHITE)




